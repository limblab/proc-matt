% plots error metrics as function of trial using trial_data struct
close all;

% get average error to each target in Baseline
utheta = unique([trial_data(get_trial_data_indices(trial_data,'epoch','bl')).target_direction]);

bl_errs = zeros(1,length(utheta));
for idx_targ = 1:length(utheta)
    idx_trial = find(get_trial_data_indices(trial_data,'epoch','bl','target_direction',utheta(idx_targ)));
    
    hand_dir = zeros(1,length(idx_trial));
    for i = 1:length(idx_trial)
        idx = idx_trial(i);
        
        hand_dir(i) = atan2( ...
        trial_data(idx).vel(trial_data(idx).idx_movement_on+15,2) - trial_data(idx).vel(trial_data(idx).idx_movement_on,2), ...
        trial_data(idx).vel(trial_data(idx).idx_movement_on+15,1) - trial_data(idx).vel(trial_data(idx).idx_movement_on,1));
    end
    
    bl_errs(idx_targ) = circular_mean(hand_dir' - utheta(idx_targ));
end

errors = zeros(1,length(trial_data));
for idx_trial = 1:length(trial_data)
    hand_dir = atan2( ...
        trial_data(idx_trial).vel(trial_data(idx_trial).idx_movement_on+15,2) - trial_data(idx_trial).vel(trial_data(idx_trial).idx_movement_on,2), ...
        trial_data(idx_trial).vel(trial_data(idx_trial).idx_movement_on+15,1) - trial_data(idx_trial).vel(trial_data(idx_trial).idx_movement_on,1));
    
    errors(idx_trial) = angleDiff(trial_data(idx_trial).target_direction,hand_dir) - bl_errs(utheta == trial_data(idx_trial).target_direction);
end
% (get_trial_data_indices(trial_data,'epoch','bl'))
p = max(bootci(1000,{@(x) circular_std(x), errors'},'alpha',0.01));

idx = abs(errors) > 1.1;

errors = errors(~get_trial_data_indices(trial_data,'epoch','bl'));
idx = idx(~get_trial_data_indices(trial_data,'epoch','bl'));


%%%
glm_encoding_model_plots;

x = errors(~idx); y = nanmean(d_vaf(~idx),1);
figure;
plot(x,y,'+')

[fit,~,~,~,s] = regress(y',[ones(size(x,2),1), x']);



% make scatter plot of r2 against errors