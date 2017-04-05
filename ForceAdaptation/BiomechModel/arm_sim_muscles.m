cf_count = 0;
for iTrial = 1:length(trial_data)
    
    if strcmpi(trial_data(iTrial).epoch,'AD')
        K = trial_data(iTrial).perturbation_info(1); % curl field constant
    else
        K = 0;
    end
    TH_c = trial_data(iTrial).perturbation_info(2); % angle of curl field application
    
    % get position/velocity and convert to meters
    p = (trial_data(iTrial).pos - repmat(pos_offset,size(trial_data(iTrial).pos,1),1) + repmat(origin_pos,size(trial_data(iTrial).pos,1),1)) / 100;
    v = trial_data(iTrial).vel / 100;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% calculate joint angles
    th = zeros(size(p,1),2);
    for t = 1:size(p,1)
        
        if left_arm
            p(t,1) = -p(t,1); % flip position
        end
        
        d = (p(t,1)^2 + p(t,2)^2 - L1^2 - L2^2)/(2*L1*L2);
        
        if d > 1
            disp('OH CRAP SOMETHING IS WEIRD IN THE POSITION');
            d = 1;
        end
        
        % there are two solutions to the quadratic, so pick one in the bound
        %th(t,2) = atan2(sqrt(1-d^2) , d);
        th(t,2) = acos(d);
        
        th(t,1) = atan2(p(t,2),p(t,1)) - atan2( (L2*sin(th(t,2))) , (L1+L2*cos(th(t,2))) );
        
        if left_arm
            th(t,1) = pi-th(t,1); % adjust angles
            th(t,2) = 2*pi-th(t,2);
            p(t,1) = -p(t,1); % flip it back
        end
    end
    
    % get joint angular velocity and acceleration
    dth = zeros(size(th,1),2);
    ddth = zeros(size(th,1),2);
    dth(:,1) = gradient(th(:,1),dt);
    dth(:,2) = gradient(th(:,2),dt);
    ddth(:,1) = gradient(dth(:,1),dt);
    ddth(:,2) = gradient(dth(:,2),dt);
    
    %%% calculate curl field force vector
    Fc = zeros(size(v,1),2);
    for t = 1:size(v,1)
        Fc(t,:) = 100 * K * [cos(TH_c) * v(t,1) + sin(TH_c)*v(t,2), -sin(TH_c)*v(t,1) + cos(TH_c) * v(t,2)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% calculate joint torques at each time bin
    
    % calculate constants for equations of motion
    R1 = 1/2*L1;
    R2 = 1/2*L2;
    I1 = 1/3*M1*(L1^2);
    I2 = 1/3*M2*(L2^2);
    
    A = I1 + I2 + M1*(R1^2) + M2*(L1^2 + R2^2);
    B = M2*L1*R2;
    C = I2 + M2*(R2^2);
    
    % loop through time and compute torques
    T = zeros(size(ddth,1),2);
    T_plan = zeros(size(ddth,1),2);
    T_force = zeros(size(ddth,1),2);
    for t = 1:size(ddth,1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % build matrix of inertial terms
        ddTerms = [A + 2*B*cos(th(t,2)), C + B*cos(th(t,2)); ...
            C + B*cos(th(t,2)),   C];
        % build matrix of coriolis terms
        dTerms = [-B*sin(th(t,2))*dth(t,2), -B*sin(th(t,2))*(dth(t,1) + dth(t,2)); ...
            B*sin(th(t,2))*dth(t,1),  0];
        
        % compute torques for this time - cross product
        % forceTorques1 = cross([p(t,:),0],[Fc(t,:),0]); % shoulder torque
        % forceTorques2 = cross([p(t,:)-[L1*cos(th(t,1)), L1*sin(th(t,1))],0],[Fc(t,:),0]); % elbow torque
        
        % jacobian method (gives same result as cross product)
        J = [-L1*sin(th(t,1)) - L2*sin(th(t,1)+th(t,2)), -L2*sin(th(t,1)+th(t,2));
            L1*cos(th(t,1)) + L2*cos(th(t,1)+th(t,2)),  L2*cos(th(t,1)+th(t,2))];
        cfTorques = J'*Fc(t,:)';
        
        T(t,:) = ddTerms * ddth(t,:)' + dTerms * dth(t,:)' + cfTorques;
        
        T_plan(t,:) = ddTerms * ddth(t,:)' + dTerms * dth(t,:)';
        T_force(t,:) = cfTorques;
    end
    
    sim_data(iTrial).kin.real.pos = p;
    sim_data(iTrial).kin.real.vel = v;
    sim_data(iTrial).torques = T;
    sim_data(iTrial).torques_plan = T_plan;
    sim_data(iTrial).torques_force = T_force;
    sim_data(iTrial).kin.real.angles = th;
    sim_data(iTrial).kin.real.dangles = dth;
    sim_data(iTrial).kin.real.ddangles = ddth;
end

%% Calculate muscle activations
disp('Calculating muscle activations...');

if simple_muscle_model
    how_many_muscles = 4;
    % This is insertion distance as proportion of segment length
    %   For now, assume muscles have same lever arm on both segments
    muscle_d = [0.02, 0.02, 0.02, 0.02]; % only for simple model
    
    for iTrial = 1:length(trial_data)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate muscle activations
        %   - scaled from 0 to 1
        %   [sh flex, sh ext, el flex, el ext]
        muscles = zeros(size(sim_data(iTrial).torques,1),4);
        
        % calculate muscle action angles as a function of time
        muscle_angles = repmat([pi/2,pi/2,pi/2,pi/2],size(sim_data(iTrial).torques,1),1);
        %muscle_angles = repmat([15,4.88,80.86,19.32]*pi/180,size(sim_data(iTrial).torques,1),1); % from Lillicrap supplementary materials
        %muscle_angles = repmat([45,30,45,30]*pi/180,size(sim_data(iTrial).torques,1),1);
        
        % Now calculate muscle force needed to cause the observed torque
        % shoulder flexion
        idx = sim_data(iTrial).torques(:,1) > 0;
        muscles(idx,1) = sim_data(iTrial).torques(idx,1)./(muscle_d(1)*sin(muscle_angles(idx,1)));
        
        % shoulder extension
        idx = sim_data(iTrial).torques(:,1) < 0;
        muscles(idx,2) = abs(sim_data(iTrial).torques(idx,1))./(muscle_d(2)*sin(muscle_angles(idx,2)));
        
        % elbow flexion
        idx = sim_data(iTrial).torques(:,2) > 0;
        muscles(idx,3) = sim_data(iTrial).torques(idx,2)./(muscle_d(3)*sin(muscle_angles(idx,3)));
        
        % elbow extension
        idx = sim_data(iTrial).torques(:,2) < 0;
        muscles(idx,4) = abs(sim_data(iTrial).torques(idx,2))./(muscle_d(4)*sin(muscle_angles(idx,4)));
        
        sim_data(iTrial).muscles = muscles;
        
    end
    
else % now it becomes an optimization problem
    how_many_muscles = 6;
    % define parameters from Lillicrap/Scott 2013
    switch how_many_muscles
        case 6
            M = [2 -2 0  0 1.5 -2   ; %moment arm matrix
                0  0 2 -2  2  -1.5];
            InsAng = 90*[1 1 1 1 1 1;
                1 1 1 1 1 1]*pi/180;
            
            PCSA = [7.6+1.9, ... % sh flex: pectoralis major, anterior deltoid
                2.0+5.2, ... % sh ext: posterior deltoid, deltoid middle
                3.7+1.7+1.5, ... % elb flex: brachialis, brachioradialis, extensor carpi radialis longus
                7.7+10.4, ... % elb ext: triceps lateral + triceps long
                4.6+2.6, ... % bi flex: biceps long, biceps short
                1.4+10.4]; % bi ext: dorsoepitrochlearis, triceps long
            
            % optimal lengths and angles
            theta0 = [15 4.88 0 0 4.5 2.12; 0 0 80.86 109.32 92.96 91.52]*pi/180;
            L0 = [7.32 3.26 6.4 4.26 5.95 4.04];
        case 4
            M = [2 -2 0  0   ; %moment arm matrix
                0  0 2 -2];
            % optimal lengths and angles
            theta0 = [15 4.88 0 0; 0 0 80.86 109.32]*pi/180;
            L0 = [7.32 3.26 6.4 4.26];
            
            PCSA = [7.6+1.9, ... % sh flex: pectoralis major, anterior deltoid
                2.0+5.2, ... % sh ext: posterior deltoid, deltoid middle
                3.7+1.7+1.5, ... % elb flex: brachialis, brachioradialis, extensor carpi radialis longus
                7.7+10.4]; % elb ext: triceps lateral + triceps long
    end
    
    params.PCSA = PCSA;
    params.InsAng = InsAng;
    
    % parameters for F-L and F-v curves
    B = 1.55;
    w = 0.81;
    Vmax = -7.39;
    cv0 = -3.21;
    cv1 = 4.17;
    bv = 0.62;
    av0 = -3.12;
    av1 = 4.21;
    av2 = -2.67;
    
    opts = optimoptions(@fmincon,'Display','off');
    
    disp('Computing muscle lengths and velocities');
    % loop along time to get muscle lengths and velocities
    for iTrial = 1:length(sim_data)
        muscle_L = zeros(length(sim_data(iTrial).torques),size(M,2));
        muscle_v = zeros(length(sim_data(iTrial).torques),size(M,2));
        for i = 1:size(M,2) % loop along muscles
            for t = 1:length(sim_data(iTrial).torques)
                muscle_L(t,i) = 1 + (M(1,i)*(theta0(1,i) - sim_data(iTrial).kin.real.angles(t,1)))/L0(i) + (M(2,i)*(theta0(2,i) - sim_data(iTrial).kin.real.angles(t,2)))/L0(i);
                muscle_v(t,i) = (M(1,i)*sim_data(iTrial).kin.real.dangles(t,1))/L0(i) + (M(2,i)*sim_data(iTrial).kin.real.dangles(t,2))/L0(i);
            end
        end
        sim_data(iTrial).muscle_L = muscle_L;
        sim_data(iTrial).muscle_v = muscle_v;
    end
    
    % now optimize muscle forces
    for iTrial = 1:length(sim_data)
        disp(['Trial ' num2str(iTrial) ' of ' num2str(length(sim_data))]);
        
        T = sim_data(iTrial).torques*100; % NOTE CONVERT TO CM BECAUSE THAT'S WHAT LILLICRAP USED
        Fv = sim_data(iTrial).muscle_v;
        Fl = sim_data(iTrial).muscle_L;
        
        muscles = zeros(length(sim_data(iTrial).torques),size(M,2));
        muscle_flags = zeros(length(sim_data(iTrial).torques),1);
        parfor t = 1:length(sim_data(iTrial).torques)
            % seed a solution for muscle forces
            u = 10*ones(1,size(M,2));
            
            if use_insertion_angles
                % use insertion angles and moment arms to get transformation relating muscle activation to torque
                
            end
            
            fuck = 1;
            switch fuck
                case 1 % no length/velocity stuff
                    A = M;
                case 2 % only length
                    A = M.*[Fl(t,:); Fl(t,:)];
                case 3 % both
                    A = M.*[Fl(t,:).*Fv(t,:); Fl(t,:).*Fv(t,:)];
            end
            
            A = A.*sin(InsAng);
            
            [muscles(t,:),~,exitflag] = fmincon(@(x) norm(x./PCSA)^2,u,[],[],A,T(t,:)',zeros(1,length(u)),Inf*ones(1,length(u)),[],opts);
            muscle_flags(t) = exitflag;
            if exitflag ~= 1
                warning('could be bad solution...')
            end
        end
        sim_data(iTrial).muscles = muscles;
        sim_data(iTrial).muscle_flags = muscle_flags;
    end
end

% First, get max torque across all trials
T_max = max(cell2mat(cellfun(@(x) max(x)',{sim_data.torques},'UniformOutput',false))',[],1);
T_min = min(cell2mat(cellfun(@(x) min(x)',{sim_data.torques},'UniformOutput',false))',[],1);
% Now, max velocity across all trials
% V_max = max(cell2mat(cellfun(@(x) max(x)',{sim_data.kin.real.vel},'UniformOutput',false))',[],1);
% V_min = min(cell2mat(cellfun(@(x) min(x)',{sim_data.kin.real.vel},'UniformOutput',false))',[],1);
% Now, max muscles across all trials
%M_max = max(cell2mat(cellfun(@(x) max(x)',{sim_data.muscles},'UniformOutput',false))',[],1);
M_max = max(prctile(cell2mat(cellfun(@(x) (x)',{sim_data.muscles},'UniformOutput',false)),[1 99],2),[],2)';
M_min = min(cell2mat(cellfun(@(x) min(x)',{sim_data.muscles},'UniformOutput',false))',[],1);

params.M_max = M_max;
params.M_min = M_min;
params.T_max = T_max;
params.T_min = T_min;
params.how_many_muscles = how_many_muscles;
