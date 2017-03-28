%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run GLM
model_names = fieldnames(cov_in);
for iFile = 1:length(session_idx)
    disp(['Fitting GLMs for File ' num2str(iFile) ' of ' num2str(length(session_idx)) '.']);
    file = session_idx(iFile);
    
    % get relevant trials
    [trial_data,params] = loadTDfiles(filepaths{iFile},func_calls{:});
    train_trials = getTDidx(trial_data,'epoch',train_info{1},'range',train_info{2});
    
    if size(trial_data(1).([in_array '_spikes']),2) == 0 || size(trial_data(1).([out_array '_spikes']),2) == 0
        error('No signals to predict or use as inputs!');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get PC projections
    if any(ismember(model_names,{'potent','null','m1_pca','pmd_pca'}))
        % get dimensionality if desired
        if isempty([in_dims, out_dims])
            disp('Getting dimensionality with Machens method.');
            switch lower(null_size)
                case 'double'
                    [~,td_temp] = getTDidx(trial_data,'epoch',{'BL','AD'});
                    td_temp = smoothSignals(td_temp,struct('signals',{{'M1_spikes','PMd_spikes'}},'kernel_SD',pn_kernel_SD));
                    td_temp = trimTD(td_temp,{'idx_movement_on',-round(30/num_bins)},{'idx_movement_on',round(50/num_bins)});
                    temp_out_dims = estimateDimensionality(td_temp,struct( ...
                        'signals',[out_array '_spikes']));
                    temp_in_dims = 2*temp_out_dims;
                otherwise
                    %trim_idx = {'idx_movement_on',-50;'idx_movement_on',50};
                    [~,td_temp] = getTDidx(trial_data,'epoch',{'BL','AD'});
                    td_temp = smoothSignals(td_temp,struct('signals',{{'M1_spikes','PMd_spikes'}},'kernel_SD',pn_kernel_SD));
                    td_temp = appendTDs( ...
                        trimTD(td_temp,{'idx_target_on',0},{'idx_target_on',round(50/num_bins)}), ...
                        trimTD(td_temp,{'idx_go_cue',0},{'idx_go_cue',round(15/num_bins)}), ...
                        trimTD(td_temp,{'idx_movement_on',0},{'idx_movement_on',round(50/num_bins)}));
                    temp_in_dims = estimateDimensionality(td_temp,struct( ...
                        'signals',[in_array '_spikes']));
                    temp_out_dims = estimateDimensionality(td_temp,struct( ...
                        'signals',[out_array '_spikes']));
            end, clear td_temp;
        else
            temp_in_dims = in_dims;
            temp_out_dims = out_dims;
        end
        
        disp('Doing some smoothing shit')
        
        td = trial_data;
        
        %trial_data = sqrtTransform(trial_data,{'M1_spikes','PMd_spikes'});
        %trial_data = smoothSignals(trial_data,struct('signals',{{'M1_spikes','PMd_spikes'}},'kernel_SD',0.1));
        
        potentNull_params = struct( ...
            'in_signals',[in_array '_spikes'], ...
            'out_signals',[out_array '_spikes'], ...
            'in_dims',temp_in_dims, ...
            'out_dims',temp_out_dims, ...
            'use_trials',{{'epoch','BL'}}, ...
            'sqrt_transform',true, ...
            'do_smoothing',true, ...
            'kernel_SD',pn_kernel_SD);
        trial_data = getPotentSpace(trial_data,potentNull_params);
        params.potentNull_params = potentNull_params;
        
        for trial = 1:length(td)
            trial_data(trial).M1_spikes = td(trial).M1_spikes;
            trial_data(trial).PMd_spikes = td(trial).PMd_spikes;
        end, clear td;
        
        trial_data = sqrtTransform(trial_data,{'PMd_spikes'});
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Do that thing where I lag kinematic signals
    if kin_lags > 0
    for trial = 1:length(trial_data)
        for i = 1:length(basic_inputs)
            temp = trial_data(trial).(basic_inputs{i});
            trial_data(trial).(basic_inputs{i}) = [temp(kin_lags+1:end,:); NaN(kin_lags,size(temp,2))];
        end
    end
    end
    
    trial_data = trimTD(trial_data,idx_start,idx_end);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit the models
    for m = 1:length(model_names)
        disp(['Calculating ' model_names{m} ' GLM...']); tic;
        [trial_data, glm_info.(model_names{m})] = getModel(trial_data,struct( ...
            'model_type','glm', ...
            'in_signals',{cov_in.(model_names{m})}, ...
            'out_signals',{{[out_array '_spikes'],'all'}}, ...
            'model_name',model_names{m}, ...
            'train_idx',train_trials, ...
            'do_lasso',do_lasso, ...
            'lasso_lambda',lasso_lambda, ...
            'lasso_alpha',lasso_alpha)); toc;
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save it
    params.glm_info = glm_info;
    save(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[filenames{iFile} '_fit.mat']),'trial_data','params')
end

