% replace relevant params struct fields with values above
params.num_neurons = num_neurons;
params.muscle_gains = muscle_gains;
params.use_muscle_model = use_muscle_model;
params.mean_lag = mean_lag;
params.std_lag = std_lag;
params.weight_distribution = weight_distribution;

%
use_models = params.use_models; if ~iscell(use_models),use_models = {use_models}; end
simple_muscle_model = params.simple_muscle_model;
how_many_muscles = params.how_many_muscles;
%         M_max = params.M_max;
M_max = max(prctile(cell2mat(cellfun(@(x) (x)',{sim_data.muscles},'UniformOutput',false)),[2.5 97.5],2),[],2)';

M_min = params.M_min;
use_muscle_model = params.use_muscle_model;

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
        
        tc_gain = [0.05 ones(1,size(sim_data(1).muscles,2))]; % OPTION 2 (ALSO CHANGE FR_GAIN)
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
                    rand_weights = 1.5*normrnd(0,1,1,size(tc_gain,2)-1);
                case 'skewgauss'
                    rand_weights = 0.25*normrnd(0.5,1,1,size(tc_gain,2)-1);
            end
            
            switch lower(use_muscle_model)
                case 'synergy' % pick flexors OR extensors
                    if rand > 0.5 % flexors
                        temp_gain = tc_gain .* [1 repmat([1 0],1,how_many_muscles/2)];
                    else % extensors
                        temp_gain = tc_gain .* [1 repmat([0 1],1,how_many_muscles/2)];
                    end
                    tc(unit,1) = temp_gain(1) .* rand(1);
                    tc(unit,2:end) = temp_gain(2:end) .* rand_weights .* muscle_gains(1:how_many_muscles);
                case 'joint'
                    if rand > 0.5 % shoulder
                        temp_gain = tc_gain .* [1 1 1 0 0 1 0];
                    else % elbow
                        temp_gain = tc_gain .* [1 0 0 1 1 0 1];
                    end
                    temp_gain = temp_gain(1:how_many_muscles+1);
                    
                    tc(unit,1) = temp_gain(1) .* rand(1);
                    tc(unit,2:end) = temp_gain(2:end) .* rand_weights .* muscle_gains(1:how_many_muscles);
                case 'single'
                    tc(unit,1) = tc_gain(1) .* rand(1);
                    m_idx = randi(how_many_muscles);
                    tc(unit,m_idx+1) = tc_gain(m_idx+1)*rand(1);
                otherwise % select from all muscles equally
                    tc(unit,1) = tc_gain(1) .* rand(1);
                    tc(unit,2:end) = tc_gain(2:end) .* rand_weights .* muscle_gains(1:how_many_muscles);
            end
            
            
            % get a random unit-specific lag
            %unit_lag(unit) = floor(normrnd(mean_lag,std_lag));
            unit_lag(unit) = randi([mean_lag-std_lag,mean_lag+std_lag],1);
            
        end
        
        unit_lag(unit_lag < 1) = 1;
        
        for iTrial = 1:length(sim_data)
            
            fr = zeros(size(sim_data(iTrial).torques,1),num_neurons);
            for unit = 1:num_neurons
                fr_gain = 1;
                n_t = size(sim_data(iTrial).torques,1);
                
                temp = 2*tc(unit,1) + fr_gain * (sum(repmat(tc(unit,2:end),n_t,1) .* ((sim_data(iTrial).muscles)./repmat(M_max,n_t,1)),2));
                
                % put cap at 1 and floor at 0
                temp(temp < tc(unit,1)) = tc(unit,1);
                temp(temp > 1) = 1;
                
                %shift neural data back by lags and pad with zeros
                temp = [temp(unit_lag(unit)+1:end); zeros(unit_lag(unit),1)];
                
                fr(:,unit) = poissrnd(temp);
            end
            sim_data(iTrial).([um '_neurons']) = fr;
        end
    else
        error('Model not recognized.');
    end
    neural_tcs.(um) = tc;
end