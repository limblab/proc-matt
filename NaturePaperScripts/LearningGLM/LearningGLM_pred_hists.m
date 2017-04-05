clc; close all; clearvars -except ANALYSIS_NAME;
dataSummary;

ANALYSIS_NAME = 'PMdM1_glm_DECENT';
load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,'start_params.mat'),'filenames');
ANALYSIS_NAME = 'PMdM1_glm_DECENT';
%filenames = filenames(cellfun(@(x) ~isempty(strfind(x,'Mihili')) | ~isempty(strfind(x,'Chewie')),filenames));

task = 'FF';
models = {'pmd'};
which_metric = 'rpr2'; %'rpr2','pr2'

bin_size = 0.05;

test_trials = {'AD',[1 10],[71,80]};
min_r2 = 0;

do_diff = false;
do_norm = false;
rem_outliers = false;


results = repmat(struct('good_cells',[],'cv',[],'pr2',[]),1,length(models));
for iModel = 1:length(models)
    fns = filenames(cellfun(@(x) ~isempty(strfind(x,['_' task '_'])),filenames));
    
    [r2] = deal([]);
    for iFile = 1:length(fns)
        tic;
        
        
        disp(['Loading File ' num2str(iFile) ' of ' num2str(length(fns)) '.']);
        load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[fns{iFile} '_fit.mat']));
        
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
    
        
    
    
    % get the values for each test trial
    [cv_ref,cv_basic] = deal([]);
    for iFile = 1:length(fns)
        load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[fns{iFile} '_cv.mat']),'cv_results');
        temp = cv_results.([models{iModel}]).(which_metric);
        if size(temp,3) > 1
            temp = squeeze(mean(temp,3));
        end
        cv_ref = cat(1,cv_ref,temp);
        
        temp = cv_results.basic.pr2;
        if size(temp,3) > 1
            temp = squeeze(mean(temp,3));
        end
        cv_basic = cat(1,cv_basic,temp);
    end
    
    % get the good cells
    good_cells = min(cv_ref,[],2) > min_r2;
    cv_ref = mean(cv_ref,2);
    cv_basic = mean(cv_basic,2);
    
    disp(['Total good cells: ' num2str(sum(good_cells)) ' / ' num2str(length(good_cells))]);
    disp(['CV pR2: ' num2str(mean(cv_ref(good_cells))) ' +/- ' num2str(std(cv_ref(good_cells)))]);
    
    
    
    

    
    if do_diff
        r2 = r2 - repmat(cv_ref',2,1);
    end
    if do_norm
        r2 = r2 ./ repmat(cv_ref',2,1);
    end
    
    r2 = r2(:,good_cells);
    
    if rem_outliers
        r2(abs(r2) > 1e5) = NaN;
        num_std = 10;
        %bad_idx = abs(r2) > repmat(num_std*std(r2,[],2),1,size(r2,2));
        bad_idx = abs(r2) > num_std*nanstd(reshape(r2,numel(r2),1));
        r2(bad_idx) = NaN;
    end
    
    x_min = min(nanmin(r2));
    x_max = max(nanmax(r2));
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