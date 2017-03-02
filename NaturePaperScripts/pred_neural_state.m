% load data
clear;
close all;
clc;

load('/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-09-15.mat');
[~,trial_data] = getTDidx(trial_data,'result','R');
trial_data = getMoveOnsetAndPeak(trial_data);


%%
[~,td] = getTDidx(trial_data,'epoch',{'BL','AD'});

td = removeBadNeurons(td,struct('min_fr',2));
td = removeBadTrials(td,struct('ranges', ...
    {{'idx_go_cue','idx_movement_on',[5,50]}}));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ESTIMATE DIMENSIONALITY
% % td_temp = appendTDs( ...
% %     truncateAndBin(td,{'idx_target_on',0},{'idx_target_on',50}), ...
% %     truncateAndBin(td,{'idx_go_cue',0},{'idx_go_cue',50}));
% td_temp = truncateAndBin(td,{'idx_go_cue',0},{'idx_go_cue',50});
% td_temp = sqrtTransform(td_temp);
% td_temp = smoothSignals(td_temp);
% % td_temp = softNormalize(td_temp);
% % td_temp = subtractConditionMean(td_temp);
% in_dims = estimateDimensionality(td_temp,struct('signal','PMd_spikes'));
% out_dims = estimateDimensionality(td_temp,struct('signal','M1_spikes'));
% clear td_temp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth and such
td = sqrtTransform(td);
td = smoothSignals(td,struct('signals',{getTDfields(td,'spikes')},'do_smoothing',true,'kernel_SD',0.1));
% td = truncateAndBin(td,{'idx_target_on',0},{'idx_trial_end',0});
td = truncateAndBin(td,{'idx_go_cue',-10},{'idx_go_cue',50});

%%
in_dims = 10;
out_dims = 5;
pca_params = struct( ...
    'in_signals','PMd_spikes', ...
    'out_signals','M1_spikes', ...
    'in_dims',in_dims, ...
    'out_dims',out_dims, ...
    'use_trials',getTDidx(td,'epoch','BL'));
[td,temp] = getPotentSpace(td,pca_params);

%%
% td = softNormalize(td);

td1 = subtractConditionMean(td(getTDidx(td,'epoch','BL')));
td2 = subtractConditionMean(td(getTDidx(td,'epoch','AD')));
td = [td1,td2];

%%
clear s;
num_dims = out_dims;
num_cv = 30;
models = {'null','potent';'null','M1';'potent','M1'};

trial_idx = getTDidx(td,'epoch','BL');
s.null_state = zeros(length(trial_idx),in_dims-out_dims);
s.potent_state = zeros(length(trial_idx),out_dims);
s.M1_state = zeros(length(trial_idx),out_dims);
s.PMd_state = zeros(length(trial_idx),out_dims+in_dims);
for i = 1:length(trial_idx)
    trial = trial_idx(i);
    idx = td(trial).idx_go_cue;
    s.null_state(i,:) = td(trial).PMdM1_null(idx,:);
    idx = td(trial).idx_movement_on;
    s.potent_state(i,:) = td(trial).PMdM1_potent(idx,:);
    s.M1_state(i,:) = td(trial).M1_pca(idx,1:out_dims);
    s.PMd_state(i,:) = td(trial).PMd_pca(idx,1:out_dims+in_dims);
end

r_cv = zeros(size(models,1),num_dims,1000);
model_b = cell(1,size(models,1));
for m = 1:size(models,1)
    invar = models{m,1};
    outvar = models{m,2};
    
    x_train_temp = s.([invar '_state']);
    y_train_temp = s.([outvar '_state']);
    
    for j = 1:1000
        % get cross validation trials
        idx = randperm(size(x_train_temp,1));
        x_cv = x_train_temp(idx(1:num_cv),:);
        y_cv = y_train_temp(idx(1:num_cv),:);
        x_train = x_train_temp(idx(num_cv+1:end),:);
        y_train = y_train_temp(idx(num_cv+1:end),:);
        
        % build model to predict potent from null
        b = zeros(size(x_train,2)+1,num_dims);
        pred_cv = zeros(size(x_cv,1),num_dims);
        for i = 1:num_dims
            b(:,i) = [ones(size(x_train,1),1), x_train]\y_train(:,i);
            pred_cv(:,i) = [ones(size(x_cv,1),1), x_cv]*b(:,i);
        end
        r_cv(m,:,j) = compute_r2(y_cv(:,1:num_dims),pred_cv);
    end
    
    b = zeros(size(x_train_temp,2)+1,num_dims);
    for i = 1:num_dims
        b(:,i) = [ones(size(x_train_temp,1),1), x_train_temp]\y_train_temp(:,i);
    end
    model_b{m} = b;
end
cv = mean(r_cv,3);

% now look at learning predictions
r = zeros(size(models,1),num_dims);

trial_idx = getTDidx(td,'epoch','AD','range',[0.5,1]);
s.null_state_ad = zeros(length(trial_idx),in_dims-out_dims);
s.potent_state_ad = zeros(length(trial_idx),out_dims);
s.M1_state_ad = zeros(length(trial_idx),out_dims);
s.PMd_state_ad = zeros(length(trial_idx),out_dims+in_dims);
for i = 1:length(trial_idx)
    trial = trial_idx(i);
    idx = td(trial).idx_go_cue;
    s.null_state_ad(i,:) = td(trial).PMdM1_null(idx,:);
    idx = td(trial).idx_movement_on;
    s.potent_state_ad(i,:) = td(trial).PMdM1_potent(idx,:);
    s.M1_state_ad(i,:) = td(trial).M1_pca(idx,1:out_dims);
    s.PMd_state_ad(i,:) = td(trial).PMd_pca(idx,1:out_dims+in_dims);
end

pred = zeros(length(trial_idx),num_dims,size(models,1));
for m = 1:size(models,1)
    b = model_b{m};
    invar = models{m,1};
    outvar = models{m,2};
    
    x_test = s.([invar '_state_ad']);
    y_test = s.([outvar '_state_ad']);
    
    for i = 1:num_dims
        pred(:,i,m) = [ones(size(x_test,1),1), x_test]*b(:,i);
    end
    r(m,:) = compute_r2(y_test(:,1:num_dims),squeeze(pred(:,:,m)));
end






% % % %%
% % % clear s;
% % % num_dims = out_dims;
% % % num_cv = 30;
% % % models = {'null','potent';'null','M1';'potent','M1'};
% % % 
% % % n_ad_trials = length(getTDidx(td,'epoch','AD'));
% % % blocks = [1,n_ad_trials];
% % % 
% % % trial_idx = getTDidx(td,'epoch','BL');
% % % s.null_state = zeros(length(trial_idx),in_dims-out_dims);
% % % s.potent_state = zeros(length(trial_idx),out_dims);
% % % s.M1_state = zeros(length(trial_idx),out_dims);
% % % s.PMd_state = zeros(length(trial_idx),out_dims+in_dims);
% % % for i = 1:length(trial_idx)
% % %     trial = trial_idx(i);
% % %     idx = td(trial).idx_go_cue;
% % %     s.null_state(i,:) = td(trial).PMdM1_null(idx,:);
% % %     idx = td(trial).idx_movement_on;
% % %     s.potent_state(i,:) = td(trial).PMdM1_potent(idx,:);
% % %     s.M1_state(i,:) = td(trial).M1_pca(idx,1:out_dims);
% % %     s.PMd_state(i,:) = td(trial).PMd_pca(idx,1:out_dims+in_dims);
% % % end
% % % 
% % % r_cv = zeros(size(models,1),num_dims,1000);
% % % for m = 1:size(models,1)
% % %     invar = models{m,1};
% % %     outvar = models{m,2};
% % %     
% % %     x_train_temp = s.([invar '_state']);
% % %     y_train_temp = s.([outvar '_state']);
% % %     
% % %     for j = 1:1000
% % %         % get cross validation trials
% % %         idx = randperm(size(x_train_temp,1));
% % %         x_cv = x_train_temp(idx(1:num_cv),:);
% % %         y_cv = y_train_temp(idx(1:num_cv),:);
% % %         x_train = x_train_temp(idx(num_cv+1:end),:);
% % %         y_train = y_train_temp(idx(num_cv+1:end),:);
% % %         
% % %         % build model to predict potent from null
% % %         b = zeros(size(x_train,2)+1,num_dims);
% % %         pred_cv = zeros(size(x_cv,1),num_dims);
% % %         for i = 1:num_dims
% % %             b(:,i) = [ones(size(x_train,1),1), x_train]\y_train(:,i);
% % %             pred_cv(:,i) = [ones(size(x_cv,1),1), x_cv]*b(:,i);
% % %         end
% % %         r_cv(m,:,j) = compute_r2(y_cv(:,1:num_dims),pred_cv);
% % %     end
% % %     
% % %     b = zeros(size(x_train_temp,2)+1,num_dims);
% % %     for i = 1:num_dims
% % %         b(:,i) = [ones(size(x_train_temp,1),1), x_train_temp]\y_train_temp(:,i);
% % %     end
% % % end
% % % cv = mean(r_cv,3);
% % % 
% % % % now look at time course
% % % r = zeros(size(models,1),num_dims,length(blocks)-1);
% % % for t = 1:length(blocks)-1
% % %     trial_idx = getTDidx(td,'epoch','AD','range',[blocks(t) blocks(t+1)]);
% % %     s.null_state_ad = zeros(length(trial_idx),in_dims-out_dims);
% % %     s.potent_state_ad = zeros(length(trial_idx),out_dims);
% % %     s.M1_state_ad = zeros(length(trial_idx),out_dims);
% % %     s.PMd_state_ad = zeros(length(trial_idx),out_dims+in_dims);
% % %     for i = 1:length(trial_idx)
% % %         trial = trial_idx(i);
% % %         idx = td(trial).idx_go_cue;
% % %         s.null_state_ad(i,:) = td(trial).PMdM1_null(idx,:);
% % %         idx = td(trial).idx_movement_on;
% % %         s.potent_state_ad(i,:) = td(trial).PMdM1_potent(idx,:);
% % %         s.M1_state_ad(i,:) = td(trial).M1_pca(idx,1:out_dims);
% % %         s.PMd_state_ad(i,:) = td(trial).PMd_pca(idx,1:out_dims+in_dims);
% % %     end
% % %     
% % %     for m = 1:size(models,1)
% % %         invar = models{m,1};
% % %         outvar = models{m,2};
% % %         
% % %         x_test = s.([invar '_state_ad']);
% % %         y_test = s.([outvar '_state_ad']);
% % %         
% % %         pred = zeros(size(x_test,1),num_dims);
% % %         for i = 1:num_dims
% % %             pred(:,i) = [ones(size(x_test,1),1), x_test]*b(:,i);
% % %         end
% % %         r(m,:,t) = compute_r2(y_test(:,1:num_dims),pred);
% % %     end
% % % end
% % % 
% % % r_norm = (r-repmat(cv,1,1,length(blocks)-1))./repmat(cv,1,1,length(blocks)-1);
% % % 
% % % close all;
% % % figure; hold all;
% % % for m = 1:size(models,1)
% % %     plot(squeeze(mean(r_norm(m,:,:),2)))
% % % end
% % % 
% % % 

