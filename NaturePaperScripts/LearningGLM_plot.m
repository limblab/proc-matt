
clc; close all; clearvars -except ANALYSIS_NAME;
dataSummary;
load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,'start_params.mat'));

models = {'pmd'}%,'potent'};

figure; hold all;
for m = 1:length(models)
    a=[];
    for iFile = 1:length(filenames)
        load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[filenames{iFile} '_eval.mat']));
        
        train_idx = params.glm_info.(models{m}).train_idx;
        ad_idx = getTDidx(trial_data,'epoch','AD','range',[1 90]);
        
%         temp = cat(2,trial_data.([models{m} '_rpr2']))';
        temp = cat(3,trial_data.([models{m} '_rpr2']));
        temp = squeeze(temp(:,1,:))';
        cv_ref = median(temp(ad_idx(end)+1:train_idx(1)-1,:),1);
        idx = cv_ref > 0.1;
        temp = temp(ad_idx,:);
        
%         temp = (temp-repmat(cv_ref,size(temp,1),1))./repmat(cv_ref,size(temp,1),1);
        
        idx = idx & all(temp > -10,1);
        
%         temp = moving_average(temp,20);
        
        a = [a, temp(:,idx)];
%                 a = [a, temp(train_idx(1)-80:train_idx(1)+50,idx)];
        
    end
    
%     temp = mean(a,2);
%     plot(temp);
    
figure;
    subplot(2,1,1);
    h1 = reshape(a(1:20,:),20*size(a,2),1);
    hist(h1,-2:0.1:1)
    set(gca,'XLim',[-2 1]);
    subplot(2,1,2);
    h2 = reshape(a(71:90,:),20*size(a,2),1);
    hist(h2,-2:0.1:1)
    set(gca,'XLim',[-2,1]);
end
