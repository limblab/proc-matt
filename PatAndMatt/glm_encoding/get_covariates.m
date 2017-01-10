function varargout = get_covariates(trial_data,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pieces together the covariates for the GLM model
% varargin:
%   params: all of them
%   start_idx: {'IDX_NAME',BINS_AFTER};
%   end_idx: same
%
% Note: params must always be an input.
%       also, start_idx must always come before end_idx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_trial = 1:length(trial_data);

if length(varargin) == 1
    params = varargin{1};
    start_idx = [];
    end_idx = [];
elseif length(varargin) == 3
    if isstruct(varargin{1})
        params = varargin{1};
        start_idx = varargin{2};
        end_idx = varargin{3};
    elseif isstruct(varargin{3})
        params = varargin{3};
        start_idx = varargin{1};
        end_idx = varargin{2};
    else
        error('Wrong inputs to get_covariates.');
    end
else
    error('Wrong inputs to get_covariates');
end

% get the relevant parameters
cov_array = params.cov_array;
pred_array = params.pred_array;
kin_signals = params.kin_signals;

if ischar(params.do_pca)
    spike_field_name = '_pca';
    spike_shift_field_name = '_pca_shift';
else
    spike_field_name = '_spikes';
    spike_shift_field_name = '_spikes_shift';
end

% get ultimate size of arrays
if isempty(start_idx)
    n_data = sum((1+cellfun(@(x) size(x,1), {trial_data(idx_trial).pos})));
else
    n_data = sum((1+[trial_data(idx_trial).(end_idx{1})]+end_idx{2}) - ([trial_data(idx_trial).(start_idx{1})]+start_idx{2}));
end

n_kin = 0;
for i = 1:length(kin_signals)
    n_kin = n_kin + size(trial_data(1).(kin_signals{i}),2);
end

n_spikes = size(trial_data(1).([cov_array spike_field_name]),2);
cov_kin = zeros(n_data,n_kin);
cov_spikes = zeros(n_data,n_spikes);

shift_data = any(~isempty(cell2mat(cellfun(@(x) strfind(x,'_shift'),fieldnames(trial_data),'Uni',0))));
if shift_data
    n_kin_shift = 0;
    for i = 1:length(kin_signals)
        n_kin_shift = n_kin_shift + size(trial_data(1).([kin_signals{i} '_shift']),2);
    end
    cov_kin_shift = zeros(n_data,n_kin_shift);
    
    if isfield(trial_data(1),spike_shift_field_name)
        n_spikes_shift = size(trial_data(1).([cov_array spike_shift_field_name]),2);
        cov_spikes_shift = zeros(n_data,n_spikes_shift);
    else
        cov_spikes_shift = [];
    end
else
    cov_kin_shift = [];
    cov_spikes_shift = [];
end

last_one = 0;
for trial = idx_trial
    if ~isempty(start_idx)
        idx = trial_data(trial).(start_idx{1})+start_idx{2}:trial_data(trial).(end_idx{1})+end_idx{2};
    else
        idx = 1:size(trial_data(trial).pos,1);
    end
    store_idx = last_one+1:last_one+length(idx);
    
    % get spiking covariate
    temp = trial_data(trial).([cov_array spike_field_name]);
    cov_spikes(store_idx,:) = temp(idx,:);
    
    % check if we want any history
    if ~isempty(cov_spikes_shift)
        temp = trial_data(trial).([cov_array spike_shift_field_name]);
        cov_spikes_shift(store_idx,:) = temp(idx,:);
    end
    
    % get kinematic covariates
    c = 0;
    for j = 1:length(kin_signals)
        temp = trial_data(trial).(kin_signals{j});
        cov_kin(store_idx,c+1:c+size(temp,2)) = temp(idx,:);
        c = c+size(temp,2);
    end
    
    % get history of kinematics if exists
    if shift_data
        c = 0;
        for j = 1:length(kin_signals)
            temp = trial_data(trial).([kin_signals{j} '_shift']);
            cov_kin_shift(store_idx,c+1:c+size(temp,2)) = temp(idx,:);
            c = c+size(temp,2);
        end
    end
    
    last_one = last_one + length(idx);
end, clear temp_kin i;

% if they aren't the same...
if ~strcmpi(cov_array,pred_array)
    
    n_spikes = size(trial_data(1).([pred_array '_spikes']),2);
    pred_spikes = zeros(n_data,n_spikes);
    
    if shift_data
        n_spikes_shift = size(trial_data(1).([pred_array '_spikes_shift']),2);
        pred_spikes_shift = zeros(n_data,n_spikes_shift);
    else
        pred_spikes_shift = [];
    end
    
    last_one = 0;
    for trial = idx_trial
        if ~isempty(start_idx)
            idx = trial_data(trial).(start_idx{1})+start_idx{2}:trial_data(trial).(end_idx{1})+end_idx{2};
        else
            idx = 1:size(trial_data(trial).pos,1);
        end
        store_idx = last_one+1:last_one+length(idx);
        
        % get spiking covariate
        temp = trial_data(trial).([pred_array '_spikes']);
        pred_spikes(store_idx,:) = temp(idx,:);
        
        % check if we want any history
        if shift_data
            temp = trial_data(trial).([pred_array '_spikes_shift']);
            pred_spikes_shift(store_idx,:) = temp(idx,:);
        end
        
        last_one = last_one + length(idx);
    end, clear i temp bin_inds fr fr_shift;
else
    pred_spikes = cov_spikes;
    pred_spikes_shift = cov_spikes_shift;
end


varargout{1} = cov_spikes;
varargout{2} = cov_spikes_shift;
varargout{3} = pred_spikes;
varargout{4} = pred_spikes_shift;
varargout{5} = cov_kin;
varargout{6} = cov_kin_shift;
