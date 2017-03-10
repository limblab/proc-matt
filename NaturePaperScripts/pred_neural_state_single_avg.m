% Predict future neural state
clear; close all; clc;
dataSummary;

in_dims  = 10;  %16;
out_dims = 5; %8;

num_cv = 10;
models = {'null','M1';'potent','M1'};

file = 70;

fname = [filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat'];
load(fullfile(rootDir,TDDir,fname),'trial_data');

[~,td] = getTDidx(trial_data,'epoch',{'BL','AD'},'result','R');

td = getMoveOnsetAndPeak(td);

td = removeBadNeurons(td,struct('min_fr',0));
td = removeBadTrials(td,struct('ranges', ...
    {{'idx_go_cue','idx_movement_on',[5,45]}}));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTIMATE DIMENSIONALITY
if isempty([in_dims,out_dims])
        td_temp = appendTDs( ...
            trimTD(td,{'idx_target_on',0},{'idx_target_on',50}), ...
            trimTD(td,{'idx_go_cue',0},{'idx_go_cue',50}));
%     td_temp = trimTD(td,{'idx_go_cue',0},{'idx_go_cue',50});
    td_temp = smoothSignals(td_temp,struct('signals',{getTDfields(td_temp,'spikes')},'sqrt_transform',true));
    td_temp = softNormalize(td_temp);
    % td_temp = subtractConditionMean(td_temp);
    in_dims = estimateDimensionality(td_temp,struct('signal','PMd_spikes'));
    out_dims = estimateDimensionality(td_temp,struct('signal','M1_spikes'));
    clear td_temp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth and such
td = smoothSignals(td,struct('signals',{getTDfields(td,'spikes')},'kernel_SD',0.1,'sqrt_transform',true));
td = trimTD(td,{'idx_target_on',0},{'idx_trial_end',-20});
% td = trimTD(td,{'idx_go_cue',-10},{'idx_movement_on',10});

td = softNormalize(td);

pca_params = struct( ...
    'in_signals','PMd_spikes', ...
    'out_signals','M1_spikes', ...
    'in_dims',in_dims, ...
    'out_dims',out_dims, ...
    'use_trials',getTDidx(td,'epoch','BL'));
[td,pca_info] = getPotentSpace(td,pca_params);

td = trimTD(td,{'idx_go_cue',-10},{'idx_go_cue',50});
td = binTD(td,10);
% td1 = subtractConditionMean(td(getTDidx(td,'epoch','BL')));
% td2 = subtractConditionMean(td(getTDidx(td,'epoch','AD')));
% td = [td1,td2];

%% BUILD MODEL
close all;
clear result;
t_lag = -1:4;
for t = 1:length(t_lag)
    t
    new_td = td;
    
    new_td = trialAverage(new_td(getTDidx(new_td,'epoch','BL','range',[0 0.9])),'target_direction');
    
    clear s;
    num_dims = out_dims;
    
    trial_idx = getTDidx(new_td,'epoch','BL');
    s.null_state = zeros(length(trial_idx),in_dims-out_dims);
    s.potent_state = zeros(length(trial_idx),out_dims);
    s.M1_state = zeros(length(trial_idx),out_dims);
    s.PMd_state = zeros(length(trial_idx),out_dims+in_dims);
    for i = 1:length(trial_idx)
        trial = trial_idx(i);
        idx = new_td(trial).idx_go_cue;
        s.null_state(i,:) = new_td(trial).PMdM1_null(idx,:);
        idx = new_td(trial).idx_go_cue+t_lag(t);
        s.potent_state(i,:) = new_td(trial).PMdM1_potent(idx,:);
        s.M1_state(i,:) = new_td(trial).M1_pca(idx,1:out_dims);
        s.PMd_state(i,:) = new_td(trial).PMd_pca(idx,1:out_dims+in_dims);
    end
    
    model_b = cell(1,size(models,1));
    for m = 1:size(models,1)
        invar = models{m,1};
        outvar = models{m,2};
        
        x_train_temp = s.([invar '_state']);
        y_train_temp = s.([outvar '_state']);
        
        x_train_temp = x_train_temp(:,1:in_dims-out_dims);
        
        b = zeros(size(x_train_temp,2)+1,num_dims);
        for i = 1:num_dims
            b(:,i) = [ones(size(x_train_temp,1),1), x_train_temp]\y_train_temp(:,i);
        end
        model_b{m} = b;
    end
    
    % CROSS VALIDATE
    new_td = td;
    new_td = trialAverage(new_td(getTDidx(new_td,'epoch','BL','range',[0.9 1])),'target_direction');
    
    r_cv = zeros(size(models,1),num_dims);
    
    trial_idx = getTDidx(new_td,'epoch','BL');
    s.null_state_cv = zeros(length(trial_idx),in_dims-out_dims);
    s.potent_state_cv = zeros(length(trial_idx),out_dims);
    s.M1_state_cv = zeros(length(trial_idx),out_dims);
    s.PMd_state_cv = zeros(length(trial_idx),out_dims+in_dims);
    for i = 1:length(trial_idx)
        trial = trial_idx(i);
        idx = new_td(trial).idx_go_cue;
        s.null_state_cv(i,:) = new_td(trial).PMdM1_null(idx,:);
        idx = new_td(trial).idx_go_cue+t_lag(t);
        s.potent_state_cv(i,:) = new_td(trial).PMdM1_potent(idx,:);
        s.M1_state_cv(i,:) = new_td(trial).M1_pca(idx,1:out_dims);
        s.PMd_state_cv(i,:) = new_td(trial).PMd_pca(idx,1:out_dims+in_dims);
    end
    
    pred_cv = zeros(length(trial_idx),num_dims,size(models,1));
    actual_cv = zeros(length(trial_idx),num_dims,size(models,1));
    for m = 1:size(models,1)
        b = model_b{m};
        invar = models{m,1};
        outvar = models{m,2};
        
        x_test = s.([invar '_state_cv']);
        y_test = s.([outvar '_state_cv']);
        
        x_test = x_test(:,1:in_dims-out_dims);
        
        actual_cv(:,:,m) = y_test(:,1:num_dims);
        for i = 1:num_dims
            pred_cv(:,i,m) = [ones(size(x_test,1),1), x_test]*b(:,i);
        end
        r_cv(m,:) = compute_r2(y_test(:,1:num_dims),squeeze(pred_cv(:,:,m)));
    end
    
    
    
    % now look at learning predictions
    % new_td = td;
    % new_td = trialAverage(new_td(getTDidx(new_td,'epoch','AD','range',[0.33 1])),'target_direction');
    new_td = [];
    blocks = 0:0.1:1;
    for i = 1:length(blocks)-1
        temp = trialAverage(td(getTDidx(td,'epoch','AD','range',blocks(i:i+1))),'target_direction');
        clear fuckthis;
        for j = 1:length(temp)
            fuckthis(j).PMdM1_null = temp(j).PMdM1_null;
            fuckthis(j).PMdM1_potent = temp(j).PMdM1_potent;
            fuckthis(j).M1_pca = temp(j).M1_pca;
            fuckthis(j).PMd_pca = temp(j).PMd_pca;
            fuckthis(j).idx_go_cue = temp(j).idx_go_cue;
            fuckthis(j).epoch = temp(j).epoch;
        end
        
        new_td = [new_td, fuckthis];
    end
    
    %     figure;
    r = zeros(size(models,1),num_dims);
    
    trial_idx = getTDidx(new_td,'epoch','AD');
    s.null_state_ad = zeros(length(trial_idx),in_dims-out_dims);
    s.potent_state_ad = zeros(length(trial_idx),out_dims);
    s.M1_state_ad = zeros(length(trial_idx),out_dims);
    s.PMd_state_ad = zeros(length(trial_idx),out_dims+in_dims);
    for i = 1:length(trial_idx)
        trial = trial_idx(i);
        idx = new_td(trial).idx_go_cue;
        s.null_state_ad(i,:) = new_td(trial).PMdM1_null(idx,:);
        idx = new_td(trial).idx_go_cue+t_lag(t);
        s.potent_state_ad(i,:) = new_td(trial).PMdM1_potent(idx,:);
        s.M1_state_ad(i,:) = new_td(trial).M1_pca(idx,1:out_dims);
        s.PMd_state_ad(i,:) = new_td(trial).PMd_pca(idx,1:out_dims+in_dims);
    end
    
    pred = zeros(length(trial_idx),num_dims,size(models,1));
    actual = zeros(length(trial_idx),num_dims,size(models,1));
    for m = 1:size(models,1)
        b = model_b{m};
        invar = models{m,1};
        outvar = models{m,2};
        
        x_test = s.([invar '_state_ad']);
        y_test = s.([outvar '_state_ad']);
        
        x_test = x_test(:,1:in_dims-out_dims);
        
        actual(:,:,m) = y_test(:,1:num_dims);
        for i = 1:num_dims
            pred(:,i,m) = [ones(size(x_test,1),1), x_test]*b(:,i);
        end
        r(m,:) = compute_r2(y_test(:,1:num_dims),squeeze(pred(:,:,m)));
        
    end
    result(t).r = r;
    result(t).r_cv = r_cv;
    result(t).pred_cv = pred_cv;
    result(t).actual_cv = actual_cv;
    result(t).pred = pred;
    result(t).actual = actual;
end

%%
plot_dims = 1:4;
plot_colors = [0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    %     0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840];

close all;
figure;
subplot(211);
imagesc(cell2mat(cellfun(@(x) x(1,:),{result.r},'uni',0)')');
subplot(212);
imagesc(cell2mat(cellfun(@(x) x(2,:),{result.r},'uni',0)')')


t =3;
actual = result(t).actual;
pred = result(t).pred;
figure;
for m = 1:size(models,1)
    for j = 1:length(plot_dims)
        subplot(length(plot_dims),size(models,1),m+(j-1)*size(models,1));
        plot(squeeze(actual(:,plot_dims(j),m)),squeeze(pred(:,plot_dims(j),m)),'o')
        title([r_cv(m,plot_dims(j)), r(m,plot_dims(j))])
    end
end