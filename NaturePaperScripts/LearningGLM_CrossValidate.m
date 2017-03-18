% Define block of trials to do CV on (probably between testing and training)

[fr_cv, fr_test, cv_basic, cv_ref] = deal([]);
for iFile = 1:length(fns)
    disp(['Running file ' num2str(iFile) ' of ' num2str(length(fns))]); tic;
    load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[fns{iFile} '_fit.mat']),'params','trial_data');
    
    train_idx = params.glm_info.(models{iModel}).train_idx;
    
    test_idx = getTDidx(trial_data,'epoch',test_trials{1},'range',test_trials{end});
    
    % get indices in this epoch that were neither test nor train
    if isempty(cv_trials)
        epoch_idx = getTDidx(trial_data,'epoch',test_trials{1});
        cv_idx = setdiff(epoch_idx,[train_idx,test_idx]);
    elseif iscell(cv_trials)
        cv_idx = getTDidx(trial_data,'epoch',cv_trials{1},'range',cv_trials{2});
    elseif ischar(cv_trials)
        epoch_idx = getTDidx(trial_data,'epoch',test_trials{1});
        cv_idx = train_idx(end)+1:epoch_idx(end);
    else
        cv_idx = max([test_idx(end)+1,train_idx(1)-cv_trials]):train_idx(1)-1;
    end
    if length(cv_idx) < 2, error('Too few CV trials!'); end
    disp(['Using ' num2str(length(cv_idx)) ' cross-validation trials.']);
    
    % Find firing rate in that block
    bin_size = trial_data(1).bin_size;
    all_spikes = cat(1,trial_data(cv_idx).(params.glm_info.(models{iModel}).out_signals{1}));
    fr_cv = [fr_cv; mean(all_spikes./bin_size,1)'];
    all_spikes = cat(1,trial_data(test_idx).(params.glm_info.(models{iModel}).out_signals{1}));
    fr_test = [fr_test; mean(all_spikes./bin_size,1)'];
    
    % do RPR2
    if isfield(params.glm_info,'basic') && strcmpi(which_metric,'rpr2')
        % Do bootstrapping of predictions in that block
        temp = evalModel(trial_data,struct( ...
            'model_type','glm', ...
            'out_signals',{params.glm_info.(models{iModel}).out_signals}, ...
            'model_name',{{'basic'}}, ...
            'eval_metric','pr2', ...
            'trial_idx',[cv_idx(1), cv_idx(end)], ...
            'num_boots',num_boots));
        cv_basic = [cv_basic; squeeze(temp)];
        
        temp = evalModel(trial_data,struct( ...
            'model_type','glm', ...
            'out_signals',{params.glm_info.(models{iModel}).out_signals}, ...
            'model_name',{{'basic',params.glm_info.(models{iModel}).model_name}}, ...
            'eval_metric','pr2', ...
            'trial_idx',[cv_idx(1), cv_idx(end)], ...
            'num_boots',num_boots));
        cv_ref = [cv_ref; squeeze(temp)];
    elseif strcmpi(which_metric,'pr2')
        temp = evalModel(trial_data,struct( ...
            'model_type','glm', ...
            'out_signals',{params.glm_info.(models{iModel}).out_signals}, ...
            'model_name',{{params.glm_info.(models{iModel}).model_name}}, ...
            'eval_metric','pr2', ...
            'trial_idx',[cv_idx(1), cv_idx(end)], ...
            'num_boots',num_boots));
        cv_ref = [cv_ref; squeeze(temp)];
    end, toc;
end


