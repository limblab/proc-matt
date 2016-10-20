function varargout = get_covariates_OG(trial_data, idx_trial, params, varargin)
% which_type: 'training', or 'testing'

if nargin > 3
    bad_cov_fr = varargin{1};
    bad_pred_fr = varargin{2};
else
    bad_cov_fr = [];
    bad_pred_fr = [];
end

% get the relevant parameters
cov_array = params.cov_array;
pred_array = params.pred_array;
start_idx = params.start_idx;
end_idx = params.end_idx;
num_samples = params.num_samples;
unit_lags = params.unit_lags;
min_fr = params.min_fr;
kin_signals = params.kin_signals;
kin_lags = params.kin_lags;

cov_fr = [];
cov_fr_shift = [];
cov_kin = [];
for i = idx_trial
    temp = full(trial_data(i).([cov_array '_spikes']));
    temp = temp(:,floor(trial_data(i).(start_idx{1}) + start_idx{2}):ceil(trial_data(i).(end_idx{1}) + end_idx{2}));
    bin_inds = 1:num_samples:size(temp,2);
    
    % fr is size bins x neurons
    fr = zeros(length(bin_inds),size(temp,1));
    for j = 1:length(bin_inds)-1
        fr(j,:) = sum(temp(:,bin_inds(j):bin_inds(j+1)),2);
    end
    
    % get kinematic covariates
    [acc,vel,pos,force] = deal([]);
    for k = 1:length(kin_lags)
        temp = trial_data(i).acc(floor(trial_data(i).(start_idx{1}) + start_idx{2})-kin_lags(k):ceil(trial_data(i).(end_idx{1}) + end_idx{2})-kin_lags(k),:);
        acc_shift = zeros(length(bin_inds),size(temp,2));
        for j = 1:length(bin_inds)-1
            acc_shift(j,:) = mean(temp(bin_inds(j):bin_inds(j+1),:),1);
        end
        acc = [acc,acc_shift];
        
        temp = trial_data(i).vel(floor(trial_data(i).(start_idx{1}) + start_idx{2})-kin_lags(k):ceil(trial_data(i).(end_idx{1}) + end_idx{2})-kin_lags(k),:);
        vel_shift = zeros(length(bin_inds),size(temp,2));
        for j = 1:length(bin_inds)-1
            vel_shift(j,:) = mean(temp(bin_inds(j):bin_inds(j+1),:),1);
        end
        vel = [vel,vel_shift];
        vel = [vel,hypot(vel_shift(:,1),vel_shift(:,2))]; % add speed
        
        temp = trial_data(i).pos(floor(trial_data(i).(start_idx{1}) + start_idx{2})-kin_lags(k):ceil(trial_data(i).(end_idx{1}) + end_idx{2})-kin_lags(k),:);
        pos_shift = zeros(length(bin_inds),size(temp,2));
        for j = 1:length(bin_inds)-1
            pos_shift(j,:) = mean(temp(bin_inds(j):bin_inds(j+1),:),1);
        end
        pos = [pos,pos_shift];
        
        temp = trial_data(i).force(floor(trial_data(i).(start_idx{1}) + start_idx{2})-kin_lags(k):ceil(trial_data(i).(end_idx{1}) + end_idx{2})-kin_lags(k),:);
        force_shift = zeros(length(bin_inds),size(temp,2));
        for j = 1:length(bin_inds)-1
            force_shift(j,:) = mean(temp(bin_inds(j):bin_inds(j+1),:),1);
        end
        force = [force,force_shift];
    end
    
    % duplicate and shift firing rates to include history as model covariates
    if unit_lags > 0
        fr_shift = zeros(size(fr,1)-unit_lags,size(fr,2),unit_lags);
        for j = 1:unit_lags
            fr_shift(:,:,j) = fr(unit_lags-(j-1):end-(j),: );
        end
        
        warning('SHOULD DO THIS SHIFT BEFORE TRUNCATING');
        fr = fr(unit_lags+1:end,:);
        acc = acc(unit_lags+1:end,:);
        vel = vel(unit_lags+1:end,:);
        pos = pos(unit_lags+1:end,:);
        force = force(unit_lags+1:end,:);
        cov_fr_shift = cat(1,cov_fr_shift,fr_shift);
    end
    
    cov_fr = [cov_fr; fr];
    
    % get kinematic covariates
    temp_kin = [];
    if any(ismember(kin_signals,'pos')), temp_kin = [temp_kin, pos]; end
    if any(ismember(kin_signals,'vel')), temp_kin = [temp_kin, vel]; end
    if any(ismember(kin_signals,'acc')), temp_kin = [temp_kin, acc]; end
    if any(ismember(kin_signals,'force')), temp_kin = [temp_kin, force]; end
    cov_kin = [cov_kin; temp_kin];
    
end, clear vel fr fr_shift j temp i bin_inds;

% check covariate firing rate matrix for low firing rates
if isempty(bad_cov_fr)
    bad_cov_fr = mean(cov_fr,1) < min_fr;
end
cov_fr(:,bad_cov_fr) = [];
if unit_lags > 0
    cov_fr_shift(:, bad_cov_fr, :) = [];
end

% if they aren't the same...
if ~strcmpi(cov_array,pred_array)
    pred_fr = [];
    pred_fr_shift = [];
    for i = idx_trial
        temp = full(trial_data(i).([pred_array '_spikes']));
        temp = temp(:,floor(trial_data(i).(start_idx{1}) + start_idx{2}):ceil(trial_data(i).(end_idx{1}) + end_idx{2}));
        bin_inds = 1:num_samples:size(temp,2);
        fr = zeros(length(bin_inds),size(temp,1));
        for j = 1:length(bin_inds)-1
            fr(j,:) = sum(temp(:,bin_inds(j):bin_inds(j+1)),2);
        end
        
        if unit_lags > 0
            fr_shift = zeros(size(fr,1)-unit_lags,size(fr,2),unit_lags);
            for j = 1:unit_lags
                fr_shift(:,:,j) = fr(unit_lags-(j-1):end-(j), : );
            end
            
            fr = fr(unit_lags+1:end,:);
            pred_fr_shift = cat(1,pred_fr_shift,fr_shift);
        end
        
        pred_fr = [pred_fr; fr];
    end, clear i temp bin_inds fr fr_shift;
    
    % check covariate firing rate matrix for low firing rates
    if isempty(bad_pred_fr)
        bad_pred_fr = mean(pred_fr,1) < min_fr;
    end
    pred_fr( :,bad_pred_fr) = [];
    pred_fr_shift(:, bad_pred_fr,: ) = [];
else
    pred_fr = cov_fr;
    pred_fr_shift = cov_fr_shift;
    bad_pred_fr = bad_cov_fr;
end


varargout{1} = cov_fr;
varargout{2} = cov_fr_shift;
varargout{3} = bad_cov_fr;
varargout{4} = pred_fr;
varargout{5} = pred_fr_shift;
varargout{6} = bad_pred_fr;
varargout{7} = cov_kin;
