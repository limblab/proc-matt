
clc; close all; clearvars -except ANALYSIS_NAME;
dataSummary;
load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,'start_params.mat'));

num_boots = 0;

model_names = fieldnames(cov_in);
if cv_blocks > 0
    for iFile = 1:length(filenames)
        tic;
        disp(['Loading File ' num2str(iFile) ' of ' num2str(length(filenames)) '.']);
        load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[filenames{iFile} '_fit.mat']));
        train_trials = getTDidx(trial_data,'epoch',train_info{1},'range',train_info{2});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cross validate the models
        
        
        % cut train_trials up into pieces
        switch lower(cv_type)
            case 'rand'
                num_cv_trials = cv_block_size;
                cv_test_trials = zeros(cv_blocks,num_cv_trials);
                cv_train_trials = zeros(cv_blocks,length(train_trials)-num_cv_trials);
                for iCV = 1:cv_blocks
                    temp = randperm(length(train_trials));
                    cv_test_trials(iCV,:) = temp(1:num_cv_trials);
                    cv_train_trials(iCV,:) = temp(num_cv_trials+1:end);
                end
            case 'fold'
                num_cv_trials = floor(length(train_trials)/cv_blocks);
                cv_test_trials = zeros(cv_blocks,num_cv_trials);
                cv_train_trials = zeros(cv_blocks,length(train_trials)-num_cv_trials);
                for iCV = 1:cv_blocks
                    cv_test_trials(iCV,:) = num_cv_trials*(iCV-1)+1:num_cv_trials*(iCV);
                    cv_train_trials(iCV,:) = setdiff(1:length(train_trials),cv_test_trials(iCV,:));
                end
        end
        
        
        td_cv = trial_data;
        num_cells = size(td_cv(1).([out_array '_spikes']),2);
        
        for m = 1:length(model_names)
            if num_boots < 2
                cv_results.(model_names{m}).pr2 = NaN(num_cells,cv_blocks);
                cv_results.(model_names{m}).rpr2 = NaN(num_cells,cv_blocks);
            else
                cv_results.(model_names{m}).pr2 = NaN(num_cells,2,cv_blocks);
                cv_results.(model_names{m}).rpr2 = NaN(num_cells,2,cv_blocks);
            end
        end
        
        % now loop along each block, fit, and test
        for iCV = 1:cv_blocks
            disp(['CV Block ' num2str(iCV) ' of ' num2str(cv_blocks)]);
            tic;
            % strip old model fits from trial_data if necessary
            for m = 1:length(model_names)
                td_cv = rmfield(td_cv,['glm_' model_names{m}]);
            end
            
            % fit all of the models on the CV trials
            for m = 1:length(model_names)
                disp(['Calculating ' model_names{m} ' GLM...']);
                [td_cv, glm_info.(model_names{m})] = getModel(td_cv,struct( ...
                    'model_type','glm', ...
                    'in_signals',{cov_in.(model_names{m})}, ...
                    'out_signals',{{[out_array '_spikes'],'all'}}, ...
                    'model_name',model_names{m}, ...
                    'train_idx',cv_train_trials(iCV,:), ...
                    'do_lasso',do_lasso, ...
                    'lasso_lambda',lasso_lambda, ...
                    'lasso_alpha',lasso_alpha));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Evaluate performance of models
            for m = 1:length(model_names)
                disp(['Evaluating ' model_names{m} ' pR2...']);
                pr2 = evalModel(td_cv,struct( ...
                    'model_type','glm', ...
                    'out_signals',{params.glm_info.(model_names{m}).out_signals}, ...
                    'model_name',params.glm_info.(model_names{m}).model_name, ...
                    'eval_metric','pr2', ...
                    'trial_idx',cv_test_trials(iCV,:), ...
                    'block_trials',true, ...
                    'num_boots',num_boots));
                if num_boots < 2
                    cv_results.(model_names{m}).pr2(:,iCV) = pr2;
                else
                    cv_results.(model_names{m}).pr2(:,:,iCV) = pr2;
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Evaluate relative performance of models
            if ismember('basic',model_names)
                temp_names = model_names(~strcmpi(model_names,'basic'));
                
                for m = 1:length(temp_names)
                    disp(['Evaluating ' temp_names{m} ' rpR2...']);
                    rpr2 = evalModel(td_cv,struct( ...
                        'model_type','glm', ...
                        'out_signals',{params.glm_info.(temp_names{m}).out_signals}, ...
                        'model_name',{{params.glm_info.basic.model_name,params.glm_info.(temp_names{m}).model_name}}, ...
                        'eval_metric','pr2', ...
                        'trial_idx',cv_test_trials(iCV,:), ...
                        'block_trials',true, ...
                        'num_boots',num_boots));
                    if num_boots < 2
                        cv_results.(temp_names{m}).rpr2(:,iCV) = rpr2;
                    else
                        cv_results.(temp_names{m}).rpr2(:,:,iCV) = rpr2;
                    end
                end, clear temp_names i m;
            end
            toc;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save it
        params.num_boots = num_boots;
        save(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[filenames{iFile} '_cv.mat']),'cv_results','params')
    end
end
