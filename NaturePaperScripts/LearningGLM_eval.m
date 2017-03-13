clc; close all; clearvars -except ANALYSIS_NAME;
dataSummary;
load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,'start_params.mat'));

model_names = fieldnames(cov_in);
for iFile = 1:length(filenames)
    tic;
    disp(['Loading File ' num2str(iFile) ' of ' num2str(length(filenames)) '.']);
    load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[filenames{iFile} '_fit.mat']));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate performance of models
    for m = 1:length(model_names)
        disp(['Evaluating ' model_names{m} ' pR2...']);
        pr2 = evalModel(trial_data,struct( ...
            'model_type','glm', ...
            'out_signals',{params.glm_info.(model_names{m}).out_signals}, ...
            'model_name',params.glm_info.(model_names{m}).model_name, ...
            'eval_metric','pr2', ...
            'trial_idx',1:length(trial_data), ...
            'num_boots',num_boots, ...
            'do_parallel',true));
        % add it to the output struct
        for trial = 1:size(pr2,1)
            trial_data(trial).([model_names{m} '_pr2']) = squeeze(pr2(trial,:,:));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate relative performance of models
    temp_names = model_names(~strcmpi(model_names,'basic'));
    for m = 1:length(temp_names)
        disp(['Evaluating ' temp_names{m} ' rpR2...']);
        rpr2 = evalModel(trial_data,struct( ...
            'model_type','glm', ...
            'out_signals',{params.glm_info.(temp_names{m}).out_signals}, ...
            'model_name',{{params.glm_info.basic.model_name,params.glm_info.(temp_names{m}).model_name}}, ...
            'eval_metric','pr2', ...
            'trial_idx',1:length(trial_data), ...
            'num_boots',num_boots, ...
            'do_parallel',true));
        % add it to the output struct
        for trial = 1:size(rpr2,1)
            trial_data(trial).([temp_names{m} '_rpr2']) = squeeze(rpr2(trial,:,:));
        end
    end, clear temp_names i m;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save it
    save(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[filenames{iFile} '_eval.mat']),'trial_data','params')
    toc;
end
