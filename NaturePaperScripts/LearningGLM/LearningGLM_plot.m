clc; close all; clearvars -except ANALYSIS_NAME;
dataSummary;

ANALYSIS_NAME = 'PMdM1_glm_test';
load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,'start_params.mat'),'filenames');

%filenames = filenames(cellfun(@(x) ~isempty(strfind(x,'Mihili')) | ~isempty(strfind(x,'Chewie')),filenames));

task = 'FF';
models = {'potent','null'};
which_metric = 'rpr2'; %'rpr2','pr2'

cv_trials = 30;
test_trials = {'AD',[1 80]};
num_boots = 100;

min_r2 = 0.0;
min_fr = 0;
r2_op = 'min'; %'mean','min'

do_norm = true;
do_diff = true;
rem_outliers = true;
mov_avg = true;

mv_avg_trials = 20;

plot_op = 'mean'; %'median','mean'

%%
plot_colors = [0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    %     0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840];
% temp_fn = filenames;
figure; hold all;
for iModel = 1:length(models)
    fns = filenames(cellfun(@(x) ~isempty(strfind(x,['_' task '_'])),filenames));
    
    % cross validate data
%     LearningGLM_CrossValidate;
    
    % get the values for each test trial
    [r2,cv_ref,cv_basic] = deal([]);
    for iFile = 1:length(fns)
        load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[fns{iFile} '_eval.mat']),'trial_data');
        load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[fns{iFile} '_cv.mat']),'cv_results');
        ad_idx = getTDidx(trial_data,'epoch',test_trials{1},'range',test_trials{2});
        r2 = [r2, cat(1,trial_data(ad_idx).([models{iModel} '_' which_metric]))];
        
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
    
    good_cells = min(cv_ref,[],2) > min_r2;
    cv_ref = mean(cv_ref,2);
    cv_basic = mean(cv_basic,2);
    % get the good cells
%     good_cells = true(size(r2,2),1);
%     fr_cells = fr_cv > min_fr & fr_test > min_fr;
%     good_cells = good_cells & fr_cells;
% 
%     switch lower(r2_op)
%         case 'min'
%             good_cells = good_cells & cv_ref(:,1) > min_r2;
% %             if strcmpi(which_metric,'rpr2')
% %                 good_cells = good_cells & cv_basic(:,1) > 0;
% %             end
%         case 'mean'
%             good_cells = good_cells & mean(cv_ref,2) > min_r2;
% %             if strcmpi(which_metric,'rpr2')
% %                 good_cells = good_cells & mean(cv_basic,2) > 0;
% %             end
%     end
    % use mean for normalizing etc
    fr_cv = 0;
    fr_cells = 0;

    
    disp(['Total good cells: ' num2str(sum(good_cells)) ' / ' num2str(length(good_cells) - sum(~fr_cells)) ' (' num2str(length(good_cells)) ')']);
    disp(['CV pR2: ' num2str(mean(cv_ref(good_cells))) ' +/- ' num2str(std(cv_ref(good_cells)))]);
    
    if do_diff
        r2 = r2 - repmat(cv_ref',size(r2,1),1);
    end
    if do_norm
        r2 = r2./repmat(cv_ref',size(r2,1),1);
    end
    
    r2 = r2(:,good_cells);
    if ~isempty(cv_basic)
        cv_basic = cv_basic(good_cells);
    end
    cv_ref = cv_ref(good_cells);
%     fr_test = fr_test(good_cells);
%     fr_cv = fr_cv(good_cells);
    
    if rem_outliers
        r2(abs(r2) > 1e5) = NaN;
        num_std = 10;
        %bad_idx = abs(r2) > repmat(num_std*std(r2,[],2),1,size(r2,2));
        bad_idx = abs(r2) > num_std*nanstd(reshape(r2,numel(r2),1));
        r2(bad_idx) = NaN;
    end
    
    if mov_avg
        r2 = moving_average(r2,mv_avg_trials);
    end
    
    m = zeros(1,size(r2,1));
    s = zeros(2,size(r2,1));
    for j = 1:size(r2,1)
        switch lower(plot_op)
            case 'mean'
                m(j) = nanmean(r2(j,:));
                s(:,j) = [m(j)-nanstd(r2(j,:))./sqrt(size(r2,2));m(j)+nanstd(r2(j,:))./sqrt(size(r2,2))];
            case 'median'
                m(j) = nanmedian(r2(j,:));
                temp = r2(j,:)';
                bs = randi(size(r2,2),size(r2,2),1000);
                s(:,j) = prctile(nanmedian(temp(bs),2),[2.5,97.5]);
        end
    end
    s = s-m(end); m = m - m(end);
    patch([1:length(m),fliplr(1:length(m))],[s(1,:),fliplr(s(2,:))],plot_colors(iModel,:),'FaceAlpha',0.3,'EdgeAlpha',0.5,'EdgeColor',plot_colors(iModel,:));
    plot(m,'LineWidth',3,'Color',plot_colors(iModel,:));
end
set(gca,'Box','off','TickDir','out','FontSize',14);