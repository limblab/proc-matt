% Predict future neural state
clear; close all; clc;
dataSummary;

num_bins = 10;

in_dims  = 16;
out_dims = 8;

num_cv = 10;
models = {'PMd','M1';'null','M1';'potent','M1'};

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
    %     td_temp = appendTDs( ...
    %         truncateAndBin(td,{'idx_target_on',0},{'idx_target_on',50}), ...
    %         truncateAndBin(td,{'idx_go_cue',0},{'idx_go_cue',50}));
    td_temp = truncateAndBin(td,{'idx_go_cue',0},{'idx_go_cue',50});
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
td = truncateAndBin(td,{'idx_target_on',0},{'idx_trial_end',0});
% td = truncateAndBin(td,num_bins,{'idx_go_cue',-10},{'idx_movement_on',10});

td = softNormalize(td);

pca_params = struct( ...
    'in_signals','PMd_spikes', ...
    'out_signals','M1_spikes', ...
    'in_dims',in_dims, ...
    'out_dims',out_dims, ...
    'use_trials',getTDidx(td,'epoch','BL'));
[td,pca_info] = getPotentSpace(td,pca_params);

td = truncateAndBin(td,{'idx_go_cue',-11},{'idx_go_cue',85});
td = truncateAndBin(td,num_bins);

%%
plot_dims = 1:3;
plot_colors = [0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    %     0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840];

close all;
figure;

num_boot = 100;
ref_idx = {'idx_go_cue',-1};
pred_idx = {'idx_peak_speed',-1};
add_offset = false;

new_td = td;

clear s;
num_dims = out_dims;

trial_idx = getTDidx(new_td,'epoch','BL');
s.null_state = zeros(length(trial_idx),in_dims-out_dims);
s.potent_state = zeros(length(trial_idx),out_dims);
s.M1_state = zeros(length(trial_idx),out_dims);
s.PMd_state = zeros(length(trial_idx),out_dims+in_dims);
for i = 1:length(trial_idx)
    trial = trial_idx(i);
    idx = new_td(trial).(ref_idx{1})+ref_idx{2};
    s.null_state(i,:) = new_td(trial).PMdM1_null(idx,:);
    idx = new_td(trial).(pred_idx{1})+pred_idx{2};
    s.potent_state(i,:) = new_td(trial).PMdM1_potent(idx,:);
    s.M1_state(i,:) = new_td(trial).M1_pca(idx,1:out_dims);
    s.PMd_state(i,:) = new_td(trial).PMd_pca(idx,1:out_dims+in_dims);
end
s.targ = [new_td(trial_idx).target_direction];

[r_cv,v_cv] = deal(zeros(size(models,1),num_dims,num_boot));
d_cv = zeros(size(models,1),num_cv,num_boot);
model_b = cell(1,size(models,1));

for m = 1:size(models,1)
    invar = models{m,1};
    outvar = models{m,2};
    
    x_train_temp = s.([invar '_state']);
    y_train_temp = s.([outvar '_state']);
    
    if ~isfield(new_td,'is_average')
        for j = 1:num_boot
            % get cross validation trials
            idx = randperm(size(x_train_temp,1));
            x_cv = x_train_temp(idx(1:num_cv),:);
            y_cv = y_train_temp(idx(1:num_cv),:);
            x_train = x_train_temp(idx(num_cv+1:end),:);
            y_train = y_train_temp(idx(num_cv+1:end),:);
            
            % build model to predict potent from null
            if add_offset
                b = zeros(size(x_train,2)+1,num_dims);
                pred_cv = zeros(size(x_cv,1),num_dims);
                for i = 1:num_dims
                    b(:,i) = [ones(size(x_train,1),1), x_train]\y_train(:,i);
                    pred_cv(:,i) = [ones(size(x_cv,1),1), x_cv]*b(:,i);
                end
            else
                b = zeros(size(x_train,2),num_dims);
                pred_cv = zeros(size(x_cv,1),num_dims);
                for i = 1:num_dims
                    b(:,i) = x_train\y_train(:,i);
                    pred_cv(:,i) =  x_cv*b(:,i);
                end
            end
            r_cv(m,:,j) = compute_r2(y_cv(:,1:num_dims),pred_cv);
            v_cv(m,:,j) = compute_vaf(y_cv(:,1:num_dims),pred_cv);
            d_cv(m,:,j) = mahal(pred_cv,y_train_temp(:,1:num_dims));
            
            subplot(2,size(models,1),m); hold all;
            for i = 1:length(plot_dims)
                plot(y_cv(:,plot_dims(i)),pred_cv(:,plot_dims(i)),'o','Color',plot_colors(i,:));
            end
            axis('tight');
            V = axis;
            plot([min(V),max(V)],[min(V),max(V)],'k--');
            axis('tight');
            set(gca,'Box','off','TickDir','out','FontSize',14);
        end
    end
    
    if add_offset
        b = zeros(size(x_train_temp,2)+1,num_dims);
        for i = 1:num_dims
            b(:,i) = [ones(size(x_train_temp,1),1), x_train_temp]\y_train_temp(:,i);
        end
    else
        b = zeros(size(x_train_temp,2),num_dims);
        for i = 1:num_dims
            b(:,i) = x_train_temp\y_train_temp(:,i);
        end
    end
    model_b{m} = b;
    
    
    title(mean(r_cv(m,plot_dims,:),3));
end
cv = mean(r_cv,3);

% now look at learning predictions
trial_idx = getTDidx(new_td,'epoch','AD','range',[0.33 1]);
s.null_state_ad = zeros(length(trial_idx),in_dims-out_dims);
s.potent_state_ad = zeros(length(trial_idx),out_dims);
s.M1_state_ad = zeros(length(trial_idx),out_dims);
s.PMd_state_ad = zeros(length(trial_idx),out_dims+in_dims);
for i = 1:length(trial_idx)
    trial = trial_idx(i);
    idx = new_td(trial).(ref_idx{1})+ref_idx{2};
    s.null_state_ad(i,:) = new_td(trial).PMdM1_null(idx,:);
    idx = new_td(trial).(pred_idx{1})+pred_idx{2};
    s.potent_state_ad(i,:) = new_td(trial).PMdM1_potent(idx,:);
    s.M1_state_ad(i,:) = new_td(trial).M1_pca(idx,1:out_dims);
    s.PMd_state_ad(i,:) = new_td(trial).PMd_pca(idx,1:out_dims+in_dims);
end
s.targ_ad = [new_td(trial_idx).target_direction];

[r,v] = deal(zeros(size(models,1),num_dims));
d = zeros(size(models,1),size(s.potent_state_ad,1));
pred = zeros(length(trial_idx),num_dims,size(models,1));
for m = 1:size(models,1)
    b = model_b{m};
    invar = models{m,1};
    outvar = models{m,2};
    
    x_test = s.([invar '_state_ad']);
    y_test = s.([outvar '_state_ad']);
    if add_offset
        for i = 1:num_dims
            pred(:,i,m) = [ones(size(x_test,1),1), x_test]*b(:,i);
        end
    else
        for i = 1:num_dims
            pred(:,i,m) = x_test*b(:,i);
        end
    end
    r(m,:) = compute_r2(y_test(:,1:num_dims),squeeze(pred(:,:,m)));
    v(m,:) = compute_vaf(y_test(:,1:num_dims),squeeze(pred(:,:,m)));
    d(m,:) = mahal(squeeze(pred(:,:,m)),y_test(:,1:num_dims));
    
    subplot(2,size(models,1),m+size(models,1)); hold all;
    for i = 1:length(plot_dims)
        plot(y_test(:,plot_dims(i)),squeeze(pred(:,plot_dims(i),m)),'o','Color',plot_colors(i,:));
    end
    title(r(m,plot_dims))
    axis('tight');
    V = axis;
    plot([min(V),max(V)],[min(V),max(V)],'k--');
    axis('tight');
    set(gca,'Box','off','TickDir','out','FontSize',14);
end
s.pred = pred;

%%
dims = 1:num_dims;
targs = unique(s.targ_ad);
% find distance of each predicted point to the cluster
d = zeros(size(s.M1_state_ad,1),size(models,1));
for m = 1:size(models,1)
    for targ = 1:length(targs)
        bl_idx = find(s.targ==targs(targ));
        ad_idx = find(s.targ_ad==targs(targ));
    
        for i = 1:length(ad_idx)
            trial = ad_idx(i);
            pred = squeeze(s.pred(trial,dims,m));
            bl = s.M1_state(bl_idx,dims);
            ad = s.M1_state_ad(ad_idx,dims);
            d(trial,m) = mahal(pred,ad);
        end
    end
end

median(d,1)

%%
dims = 1:num_dims;
figure;
[r,v] = deal(zeros(1,size(models,1)));
for m = 1:size(models,1)
    subplot(1,size(models,1),m); hold all;
    
    pred = squeeze(s.pred(:,dims,m));
    actual = s.M1_state_ad(:,dims);
    
    pred = pred./repmat(range(actual,1),size(pred,1),1);
    actual = actual./repmat(range(actual,1),size(pred,1),1);
    
    
    pred = reshape(pred,numel(pred),1);
    actual = reshape(actual,numel(actual),1);
    
    
    plot(actual,pred,'o');
    axis('tight');
    V = axis;
    plot([min(V),max(V)],[min(V),max(V)],'k--');
    axis('tight');
    set(gca,'Box','off','TickDir','out','FontSize',14);
    r(m) = compute_r2(actual,pred);
    v(m) = compute_vaf(actual,pred);
    title([r(m),mean(reshape(squeeze(mean(r_cv(m,dims,:),3)),numel(r_cv(m,dims,1)),1))]);
end
