% Predict future neural state
clear; close all; clc;
dataSummary;

num_bins = 5;

in_dims  = 8;  %16;
out_dims = 4; %8;

num_cv = 30;
models = {'null','potent';'null','M1';'potent','M1'};

tasks = {'CO'};
perts = {'FF'};

session_idx = find( ...
    ismember(filedb.Task,tasks) & ...
    ismember(filedb.Perturbation,perts) & ...
    ~cellfun(@isempty,filedb.FileNames) & ...
    ~ismember(filedb.Monkey,'MrT') & ...
    ~(ismember(filedb.Monkey,'Mihili') & datenum(filedb.Date) > datenum('2015-01-01')) & ...
    cellfun(@(x) all(ismember({'M1','PMd'},x)),filedb.Arrays));

for iFile = 1:length(session_idx)
    file = session_idx(iFile);
    disp(['File ' num2str(iFile) ' of ' num2str(length(session_idx))]);
    
    fname = [filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat'];
    load(fullfile(rootDir,TDDir,fname),'trial_data');
    if strcmpi(filedb.Task{file},'rt')
        trial_data = getRWMovements(trial_data);
    end
    
    [~,td] = getTDidx(trial_data,'epoch',{'BL','AD'},'result','R');
    
    td = getMoveOnsetAndPeak(td);
    
    td = removeBadNeurons(td,struct('min_fr',1));
    td = removeBadTrials(td,struct('ranges', ...
        {{'idx_go_cue','idx_movement_on',[5,45]}}));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ESTIMATE DIMENSIONALITY
    if isempty([in_dims,out_dims])
        td_temp = appendTDs( ...
            truncateAndBin(td,{'idx_target_on',0},{'idx_target_on',50}), ...
            truncateAndBin(td,{'idx_go_cue',0},{'idx_go_cue',50}));
        %         td_temp = truncateAndBin(td,{'idx_go_cue',0},{'idx_go_cue',50});
        td_temp = sqrtTransform(td_temp);
        td_temp = smoothSignals(td_temp,struct('signals',{getTDfields(td_temp,'spikes')}));
        % td_temp = softNormalize(td_temp);
        % td_temp = subtractConditionMean(td_temp);
        in_dims = estimateDimensionality(td_temp,struct('signal','PMd_spikes'));
        out_dims = estimateDimensionality(td_temp,struct('signal','M1_spikes'));
        clear td_temp
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % smooth and such
    td = smoothSignals(td,struct('signals',{getTDfields(td,'spikes')},'kernel_SD',0.1,'sqrt_transform',true));
    % td = truncateAndBin(td,{'idx_target_on',0},{'idx_trial_end',0});
    td = truncateAndBin(td,num_bins,{'idx_go_cue',-50},{'idx_go_cue',60});
    
    pca_params = struct( ...
        'in_signals','PMd_spikes', ...
        'out_signals','M1_spikes', ...
        'in_dims',in_dims, ...
        'out_dims',out_dims, ...
        'use_trials',getTDidx(td,'epoch','BL'));
    [td,pca_info] = getPotentSpace(td,pca_params);
    
    % td = softNormalize(td);
    %         td1 = subtractConditionMean(td(getTDidx(td,'epoch','BL')));
    %         td2 = subtractConditionMean(td(getTDidx(td,'epoch','AD')));
    %         td = [td1,td2];
    
    clear s;
    num_dims = out_dims;
    
    trial_idx = getTDidx(td,'epoch','BL');
    s.null_state = zeros(length(trial_idx),in_dims-out_dims);
    s.potent_state = zeros(length(trial_idx),out_dims);
    s.M1_state = zeros(length(trial_idx),out_dims);
    s.PMd_state = zeros(length(trial_idx),out_dims+in_dims);
    for i = 1:length(trial_idx)
        trial = trial_idx(i);
        idx = td(trial).idx_go_cue-ceil(5/num_bins);
        s.null_state(i,:) = td(trial).PMdM1_null(idx,:);
        idx = td(trial).idx_movement_on;
        s.potent_state(i,:) = td(trial).PMdM1_potent(idx,:);
        s.M1_state(i,:) = td(trial).M1_pca(idx,1:out_dims);
        s.PMd_state(i,:) = td(trial).PMd_pca(idx,1:out_dims+in_dims);
    end
    
    r_cv = zeros(size(models,1),num_dims,1000);
    model_b = cell(1,size(models,1));
    clear err_cv;
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
            % compute distance from truth
            for i = 1:size(pred_cv,1)
                err_cv(m,i,j) = norm(pred_cv(i,:)-y_cv(i,1:num_dims));
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
    
    trial_idx = getTDidx(td,'epoch','AD','range',[0 1]);
    s.null_state_ad = zeros(length(trial_idx),in_dims-out_dims);
    s.potent_state_ad = zeros(length(trial_idx),out_dims);
    s.M1_state_ad = zeros(length(trial_idx),out_dims);
    s.PMd_state_ad = zeros(length(trial_idx),out_dims+in_dims);
    for i = 1:length(trial_idx)
        trial = trial_idx(i);
        idx = td(trial).idx_go_cue-ceil(5/num_bins);
        s.null_state_ad(i,:) = td(trial).PMdM1_null(idx,:);
        idx = td(trial).idx_movement_on;
        s.potent_state_ad(i,:) = td(trial).PMdM1_potent(idx,:);
        s.M1_state_ad(i,:) = td(trial).M1_pca(idx,1:out_dims);
        s.PMd_state_ad(i,:) = td(trial).PMd_pca(idx,1:out_dims+in_dims);
    end
    
    pred = zeros(length(trial_idx),num_dims,size(models,1));
    actual = zeros(length(trial_idx),num_dims,size(models,1));
    for m = 1:size(models,1)
        b = model_b{m};
        invar = models{m,1};
        outvar = models{m,2};
        
        x_test = s.([invar '_state_ad']);
        y_test = s.([outvar '_state_ad']);
        actual(:,:,m) = y_test(:,1:num_dims);
        for i = 1:num_dims
            pred(:,i,m) = [ones(size(x_test,1),1), x_test]*b(:,i);
        end
        r(m,:) = compute_r2(y_test(:,1:num_dims),squeeze(pred(:,:,m)));
    end
    
    results(iFile).cv = cv;
    results(iFile).r = r;
    results(iFile).b = b;
    results(iFile).pred = pred;
    results(iFile).actual = actual;
    results(iFile).pca_info = pca_info;
    results(iFile).err_cv = err_cv;
end

%%
all_d = []; all_d_cv = [];
for i = 1:length(session_idx)
    d = zeros(size(models,1),size(results(i).pred,1));
    for m = 1:size(models,1)
        for k = 1:size(results(i).pred,1)
            d(m,k) = norm(squeeze(results(i).pred(k,:,m)) - squeeze(results(i).actual(k,:,m)));
        end
    end
    all_d = cat(2,all_d,d);
    all_d_cv = cat(2,all_d_cv,mean(results(i).err_cv,3));
end
