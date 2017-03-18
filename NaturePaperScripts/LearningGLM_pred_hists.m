ANALYSIS_NAME = 'PMdM1_glm_FF';

clc; close all; clearvars -except ANALYSIS_NAME;
dataSummary;
load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,'start_params.mat'));

filenames = filenames(cellfun(@(x) ~isempty(strfind(x,'_FF_')),filenames));

which_metric = 'rpr2';

test_trials = {'AD',[1 10],[71,80]};
cv_trials = {'AD',[81 90]};
min_r2 = 0.01;
min_fr = 3;
r2_op = 'mean';

do_diff = true;
do_norm = true;
rem_outliers = true;
models = {'pmd'};

num_boots = 100;

results = repmat(struct('good_cells',[],'cv',[],'pr2',[]),1,length(models));
for iModel = 1:length(models)
    
    LearningGLM_CrossValidate;
    
    [r2] = deal([]);
    for iFile = 1:length(filenames)
        tic;
        
        
        disp(['Loading File ' num2str(iFile) ' of ' num2str(length(filenames)) '.']);
        load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[filenames{iFile} '_fit.mat']));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate performance of models
        switch lower(which_metric)
            case 'pr2'
                model_ins = {params.glm_info.(models{iModel}).model_name};
            case 'rpr2'
                model_ins = {'basic',params.glm_info.(models{iModel}).model_name};
        end
        
        test_idx = getTDidx(trial_data,'epoch',test_trials{1},'range',test_trials{2});
        r2_e = evalModel(trial_data,struct( ...
            'model_type','glm', ...
            'out_signals',{params.glm_info.(models{iModel}).out_signals}, ...
            'model_name',{model_ins}, ...
            'eval_metric','pr2', ...
            'trial_idx',[test_idx(1),test_idx(end)], ...
            'num_boots',0));
        test_idx = getTDidx(trial_data,'epoch',test_trials{1},'range',test_trials{3});
        r2_l = evalModel(trial_data,struct( ...
            'model_type','glm', ...
            'out_signals',{params.glm_info.(models{iModel}).out_signals}, ...
            'model_name',{model_ins}, ...
            'eval_metric','pr2', ...
            'trial_idx',[test_idx(1),test_idx(end)], ...
            'num_boots',0));
        r2 = [r2, cat(1,r2_e,r2_l)];
    end
    
    % get the good cells
    good_cells = true(size(r2,2),1);
    fr_cells = fr_cv > min_fr & fr_test > min_fr;
    good_cells = good_cells & fr_cells;
    
    switch lower(r2_op)
        case 'min'
            cv_ref = cv_ref(:,1);
            if strcmpi(which_metric,'rpr2')
                cv_basic = cv_basic(:,1);
            end
        case 'mean'
            cv_ref = mean(cv_ref,2);
            if strcmpi(which_metric,'rpr2')
                cv_basic = mean(cv_basic,2);
            end
    end
    good_cells = good_cells & cv_ref > min_r2;
    if strcmpi(which_metric,'rpr2')
        good_cells = good_cells & cv_basic > min_r2;
    end
    
    disp(['Total good cells: ' num2str(sum(good_cells)) ' / ' num2str(length(good_cells) - sum(~fr_cells)) ' (' num2str(length(good_cells)) ')']);
    
    if do_diff
        r2 = r2 - repmat(cv_ref',2,1);
    end
    if do_norm
        r2 = r2 ./ repmat(cv_ref',2,1);
    end
    
    r2 = r2(:,good_cells);
    
    if rem_outliers
        r2(abs(r2) > 1e5) = NaN;
        num_std = 3;
        %bad_idx = abs(r2) > repmat(num_std*std(r2,[],2),1,size(r2,2));
        bad_idx = abs(r2) > num_std*nanstd(reshape(r2,numel(r2),1));
        r2(bad_idx) = NaN;
    end
    
    x_min = min(nanmin(r2));
    x_max = max(nanmax(r2));
    bin_size = 0.02;
    bins = x_min:bin_size:x_max;
    
    [~,p] = ttest2(r2(1,:),r2(2,:));
    figure; hold all;
    hist(r2(1,:),bins);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r');
    set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[x_min x_max]);
    hist(r2(2,:),bins);
    h = findobj(gca,'Type','patch');
    set(h,'EdgeColor','w','FaceAlpha',0.7,'EdgeAlpha',0.7);
    set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[x_min x_max]);
    title(p);
    
end

% [~,p] = ttest2(results(2).pr2(1,results(2).good_cells),results(3).pr2(1,results(3).good_cells))