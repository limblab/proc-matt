close all;
clear;
clc;

filenames = { ...
                'Chewie_CO_FF_2013-10-22', ...
                'Chewie_CO_FF_2013-10-23', ...
                'Chewie_CO_FF_2013-10-31', ...
                'Chewie_CO_FF_2013-11-01', ...
                'Chewie_CO_FF_2013-12-03', ...
                'Chewie_CO_FF_2013-12-04', ...
                'Chewie_CO_FF_2015-07-01', ...
                'Chewie_CO_FF_2015-07-03', ...
                'Mihili_CO_FF_2014-02-03', ...
                'Mihili_CO_FF_2014-02-17', ...
                'Mihili_CO_FF_2014-02-18', ...
                'Mihili_CO_FF_2014-03-07', ...
                'Mihili_CO_FF_2015-06-11', ...
                'Mihili_CO_FF_2015-06-17', ...
                'Chewie_CO_FF_2015-06-29', ...
                'Chewie_CO_FF_2015-06-30', ...
                'Chewie_CO_FF_2015-07-06', ...
                'Chewie_CO_FF_2015-07-07', ...
                'Chewie_CO_FF_2015-07-08', ...
                'Mihili_CO_FF_2015-06-10', ...
                'Mihili_CO_FF_2015-06-15', ...
    'Mihili_CO_FF_2015-06-16', ...
    };

for iFile = 1:length(filenames)
    
    clearvars -except filenames iFile;
    close all;
    clc;
    iFile
    
    root_dir = 'F:\trial_data_files\';
    
    %% load in trial data
    load([root_dir filenames{iFile} '.mat']);
    
    %%
    use_models = {'muscle'};
    do_poisson = true;
    left_arm = true; % flips from right arm to left arm
    
    %%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%
    pos_offset = [1,-31]; % position offset from behavior
    origin_pos = [0 20];
    dt = 0.01;
    
    %%% biomech model
    % use Cheng/Scott 2000 for parameters
    switch lower(trial_data(1).monkey)
        case 'chewie'
            body_weight = 13; % body mass in kg
            M1 = 34.4*body_weight/1000; % from Scott 2000
            M2 = 25.2*body_weight/1000;
            L1 = 0.19;
            L2 = 0.22;
            %             M1 = 0.7; % mass of upper arm in kg
            %             M2 = 0.7; % mass of lower arm in kg
            %             L1 = 0.17; % length of upper arm in m
            %             L2 = 0.17; % length of lower arm in m
        case 'mihili'
            body_weight = 9;
            M1 = 34.4*body_weight/1000;
            M2 = 25.2*body_weight/1000;
            L1 = 0.17;
            L2 = 0.20;
            %             M1 = 0.55; % mass of upper arm in kg
            %             M2 = 0.55; % mass of lower arm in kg
            %             L1 = 0.14; % length of upper arm in m
            %             L2 = 0.14; % length of lower arm in m
    end
    
    TH_1_min = 0;  % minimum shoulder angle
    TH_1_max = pi; % maximum shoulder angle
    TH_2_min = 0;  % minimum elbow angle
    TH_2_max = pi; % maximum elbow angle
    
    %%% muscle model
    % This is insertion distance as proportion of segment length
    %   For now, assume muscles have same lever arm on both segments
    muscle_d = [0.02, 0.02, 0.02, 0.02]; % 2cm insertion distance
    
    %%% neural activity model
    num_neurons = 500;
    muscle_gains = [1,1,1,1]; % gain terms for [sh flex, sh ext, el flex, el ext]
    use_muscle_model = 'synergy'; %'single','joint','synergy','all'
    weight_distribution = 'gaussian';
    mean_lag = 10; % mean lag in bins, that spikes will be shifted before tuning
    std_lag = 3; % how many bins of standard deviation
    
    if ~iscell(use_models)
        use_models = {use_models};
    end
    
    % build params struct
    params = struct('use_models',use_models,'do_poisson',do_poisson,'pos_offset',pos_offset,'origin_pos',origin_pos,'dt',dt,'M1',M1,'M2',M2,'L1',L1,'L2',L2, ...
        'muscle_d',muscle_d,'num_neurons',num_neurons,'muscle_gains',muscle_gains,'use_muscle_model',use_muscle_model,'mean_lag',mean_lag,'std_lag',std_lag);
    
    %% get behavioral data
    %     epochs = unique({trial_data.epoch});
    %
    %     idx = find(strcmpi({trial_data.epoch},'AD'));
    %
    %     trial_err = zeros(length(idx),1);
    %
    %     v = trial_data(idx(1)).vel;
    %     move_idx = trial_data(idx(1)).idx_movement_on:trial_data(idx(1)).idx_movement_on+10;
    %     trial_err(1) = angleDiff(trial_data(idx(1)).target_direction, atan2(v( move_idx(end),2)-v( move_idx(1),2),v( move_idx(end),1)-v( move_idx(1),1)), true, true);
    %
    %     for iTrial = 2:length(idx)
    %         v = trial_data(idx(iTrial)).vel;
    %         move_idx = trial_data(idx(iTrial)).idx_movement_on:trial_data(idx(iTrial)).idx_movement_on+10;
    %         trial_err(iTrial) = angleDiff(trial_data(idx(iTrial)).target_direction, atan2(v(move_idx(end),2)-v(move_idx(1),2),v(move_idx(end),1)-v(move_idx(1),1)), true, true);
    %     end
    %     % fit an exponential
    %     f = fit((1:length(idx))',trial_err,'exp1');
    %     cf_errs = f.a*exp(f.b*(1:length(idx)));
    %
    %     cf_errs_diff = [0, cumsum(diff(cf_errs))];
    %
    %     cf_direction = sign(trial_data(idx(1)).perturbation_info(2));
    %
    
    %% loop along trials and do modeling
    disp('Calculating kinematics and dynamics...');
    
    cf_count = 0;
    for iTrial = 1:length(trial_data)
        
        if strcmpi(trial_data(iTrial).epoch,'AD')
            K = trial_data(iTrial).perturbation_info(1); % curl field constant
        else
            K = 0;
        end
        TH_c = trial_data(iTrial).perturbation_info(2); % angle of curl field application
        
        % get position and convert to meters
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
            Fc(t,:) = 100 * K * [cos(TH_c)*v(t,1) + sin(TH_c)*v(t,2), sin(TH_c)*v(t,1) - cos(TH_c) * v(t,2)];
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
            if 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % build matrix of inertial terms
                ddTerms = [A + 2*B*cos(th(t,2)), C + B*cos(th(t,2)); ...
                    C + B*cos(th(t,2)),   C];
                % build matrix of coriolis terms
                dTerms = [-B*sin(th(t,2))*dth(t,2), -B*sin(th(t,2))*(dth(t,1) + dth(t,2)); ...
                    B*sin(th(t,2))*dth(t,1),  0];
                
                % compute torques for this time
                %T(t,:) = ddTerms * ddth(t,:)' + dTerms * dth(t,:)' - fxTerms*Fc(t,1) + fyTerms*Fc(t,2);
                forceTorques1 = cross([p(t,:),0],[Fc(t,:),0]); % shoulder torque
                forceTorques2 = cross([p(t,:)-[L1*cos(th(t,1)), L1*sin(th(t,1))],0],[Fc(t,:),0]); % elbow torque
                
                T(t,:) = ddTerms * ddth(t,:)' + dTerms * dth(t,:)' + [forceTorques1(3);forceTorques2(3)];
                
                T_plan(t,:) = ddTerms * ddth(t,:)' + dTerms * dth(t,:)';
                T_force(t,:) = [forceTorques1(3),forceTorques2(3)];
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % THIS METHOD ASSUMES SAME MASS AND LENGTH
                % % % %                 T(t,1) = (1/2*M1*L1^2) * ( ddth(t,1)*2*(5/3 + cos(th(t,2))) + ddth(t,2)*(2/3 + cos(th(t,2))) + sin(th(t,2))*dth(t,2)*(2*dth(t,1)+dth(t,2)) );
                % % % %                 T(t,2) = (1/2*M2*L2^2) * ( ddth(t,1)*(2/3 + cos(th(t,2))) + ddth(t,2)*2/3 + sin(th(t,2))*(dth(t,1)^2) );
                % % % %
                % % % %                 T_force(t,:) = -fxTerms*Fc(t,1) + fyTerms*Fc(t,2);
                % % % %
                % % % %                 T_plan(t,:) = T(t,:);
                % % % %                 T(t,:) = T(t,:) + T_force(t,:);
            end
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
    
    % First, get max torque across all trials
    T_max = max(cell2mat(cellfun(@(x) max(x)',{sim_data.torques},'UniformOutput',false))',[],1);
    T_min = min(cell2mat(cellfun(@(x) min(x)',{sim_data.torques},'UniformOutput',false))',[],1);
    % Now, max velocity across all trials
    % V_max = max(cell2mat(cellfun(@(x) max(x)',{sim_data.kin.real.vel},'UniformOutput',false))',[],1);
    % V_min = min(cell2mat(cellfun(@(x) min(x)',{sim_data.kin.real.vel},'UniformOutput',false))',[],1);
    % Now, max muscles across all trials
    M_max = max(cell2mat(cellfun(@(x) max(x)',{sim_data.muscles},'UniformOutput',false))',[],1);
    M_min = min(cell2mat(cellfun(@(x) min(x)',{sim_data.muscles},'UniformOutput',false))',[],1);
    
    %% Now, generate neural activity for each trial
    
    for iModel = 1:length(use_models)
        um = use_models{iModel};
        
        if strcmpi(um(1:3),'kin')
            error('kin is outdated. needs updating.');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % generate neurons from torque model
        elseif strcmpi(um,'torque')
            error('torque is outdated. needs updating');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % generate neurons from muscle model
        elseif strcmpi(um,'muscle')
            disp('Generating firing rates from muscle activity...');
            
            tc_gain = [0.05 1 1 1 1]; % OPTION 2 (ALSO CHANGE FR_GAIN)
            m_w_max = 1;
            m_w_min = -1;
            
            % linear combination of all muscles
            tc = zeros(num_neurons,size(tc_gain,2));
            unit_lag = zeros(1,num_neurons);
            for unit = 1:num_neurons
                switch lower(weight_distribution)
                    case 'uniform'
                        rand_weights = (m_w_min + (m_w_max-m_w_min).*rand(1,size(tc_gain,2)-1));
                    case 'gaussian'
                        rand_weights = normrnd(0,1,1,size(tc_gain,2)-1);
                end
                
                switch lower(use_muscle_model)
                    case 'synergy' % pick flexors OR extensors
                        if rand > 0.5 % flexors
                            temp_gain = tc_gain .* [1 1 0 1 0];
                        else % extensors
                            temp_gain = tc_gain .* [1 0 1 0 1];
                        end
                        tc(unit,1) = temp_gain(1) .* rand(1);
                        tc(unit,2:end) = temp_gain(2:end) .* rand_weights;
                    case 'joint'
                        if rand > 0.5 % shoulder
                            temp_gain = tc_gain .* [1 1 1 0 0];
                        else % elbow
                            temp_gain = tc_gain .* [1 0 0 1 1];
                        end
                        tc(unit,1) = temp_gain(1) .* rand(1);
                        tc(unit,2:end) = temp_gain(2:end) .* rand_weights;
                    case 'single'
                        tc(unit,1) = tc_gain(1) .* rand(1);
                        m_idx = randi(4);
                        tc(unit,m_idx+1) = tc_gain(m_idx+1)*rand(1);
                    otherwise % select from all muscles equally
                        tc(unit,1) = tc_gain(1) .* rand(1);
                        tc(unit,2:end) = tc_gain(2:end) .* rand_weights;
                end
                
                % get a random unit-specific lag
                %unit_lag(unit) = floor(normrnd(mean_lag,std_lag));
                unit_lag(unit) = randi([mean_lag-std_lag,mean_lag+std_lag],1);
                
            end
            
            unit_lag(unit_lag < 1) = 1;
            
            for iTrial = 1:length(sim_data)
                fr = zeros(size(sim_data(iTrial).torques,1),num_neurons);
                for unit = 1:num_neurons
                    if do_poisson
                        fr_gain = 1;
                        
                        temp = 2*tc(unit,1) + fr_gain * (sum(repmat(tc(unit,2:end),size(sim_data(iTrial).torques,1),1) .* ((sim_data(iTrial).muscles)./repmat(M_max,size(sim_data(iTrial).torques,1),1)),2));
                        % put cap at 1 and floor at 0
                        temp(temp < tc(unit,1)) = tc(unit,1);
                        temp(temp > 1) = 1;
                        
                        %shift neural data back by lags and pad with zeros
                        temp = [temp(unit_lag(unit)+1:end); zeros(unit_lag(unit),1)];
                        fr(:,unit) = poissrnd(temp);
                    else
                        fr_gain = 100;
                        temp = fr_gain * (tc(unit,1) + sum(repmat(tc(unit,2:end),size(sim_data(iTrial).torques,1),1) .* (sim_data(iTrial).muscles./repmat(M_max,size(sim_data(iTrial).torques,1),1)),2));
                        temp(temp < 0) = 0;
                        
                        %shift neural data back by lags and pad with zeros
                        temp = [temp(unit_lag(unit)+1:end); zeros(unit_lag(unit),1)];
                        fr(:,unit) = temp;
                    end
                end
                sim_data(iTrial).([um '_neurons']) = fr;
            end
        else
            error('Model not recognized.');
        end
        neural_tcs.(um) = tc;
    end
    
    %%
    arm_sim_tuning;
    %
    %     save([root_dir filenames{iFile} '_results.mat'],'sw_data','tc_data','neural_tcs','params');
    
end

%% Plot the simulation
if 1
    %     close all;
    params.dt = dt;
    params.M1 = M1; % mass of upper arm in kg
    params.M2 = M2; % mass of lower arm in kg
    params.L1 = L1; % length of upper arm in m
    params.L2 = L2; % length of lower arm in m
    params.pos_offset = pos_offset; % position offset from behavior
    params.origin_pos = origin_pos;
    params.dt = dt;
    params.M_max = M_max;
    params.M_min = M_min;
    params.T_max = T_max;
    params.T_min = T_min;
    params.num_neurons = num_neurons;
    params.comp_blocks = [1 2 3];
    
%     params.type = 'time_signals';
        params.type = 'freeze_video';
    params.resolution = 5;
    params.signals = {'vel','torques','muscle_neurons'};
    params.kin_model = 'real';
    
    params.torques_lim = [-0.5,0.25];
    params.vel_lim = [-0.3,0.1];
    
% I=cell2mat({trial_data.target_direction})==trial_data(59).target_direction & strcmpi({trial_data.epoch},'BL'); idx = find(I);
    %         idx = 43; % baseline
    %     idx = 376; % decent ad in same direction
    
    %     idx = [56,402];    
    idx = [59,438];
    
    for iTrial = idx
        iTrial
        params.idx_start = trial_data(iTrial).idx_go_cue;
        params.idx_end = trial_data(iTrial).idx_reward;
        plot_arm_sim(sim_data(iTrial),params,tc_data);
        if length(idx) > 4
            pause;
            close all;
        end
    end
end
