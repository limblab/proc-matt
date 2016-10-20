function [y, x_full, x_basic] = glm_prep_inputs(trial_data,unit,trial_idx,varargin)
% trial_data: obvious
% unit: index of predicted unit in pred_array
% trial_idx: which trials in trial_data to include
% varargin:
%   params: all of them
%   start_idx: {'IDX_NAME',BINS_AFTER};
%   end_idx: same
%
% Note: params must always be an input.
%       also, start_idx must always come before end_idx
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

cov_array = params.cov_array;
pred_array = params.pred_array;
do_pca = params.do_pca;
do_rcb = params.do_rcb;
do_all_history = params.do_all_history;
unit_lags = params.unit_lags;
kin_lags = params.kin_lags;
kin_signals = params.kin_signals;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get covariates
[cov_spikes, cov_spikes_shift, pred_spikes, pred_spikes_shift, cov_kin, cov_kin_shift] = get_covariates(trial_data,trial_idx,start_idx,end_idx,params);

[x_self, x_kin] = deal([]);
idx = true(1,size(cov_spikes,2));
if strcmpi(cov_array,pred_array)
    % take out this neuron and compute the new PCA
    idx(unit) = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% piece together covariates
y = pred_spikes(:,unit);


if do_pca % get PCA space
    % square root transform 4 lyfe
    x_spikes = sqrt(cov_spikes(:,idx))*params.pca_w;
    if ~ischar(params.pca_dims)
        x_spikes = x_spikes(:,params.pca_dims);
    end
else % just use spike counts
    x_spikes = cov_spikes(:,idx);
end


% add shifted data as covariates if desired
if unit_lags > 0
    % add time shifted other neurons (NO PCA AT THE MOMENT)
    if do_all_history && ~do_pca
        idx_shift = repmat(idx,1,unit_lags);
        x_spikes = [x_spikes, cov_spikes_shift(:,idx_shift)];
    end
    x_self = zeros(size(pred_spikes,1),unit_lags);
    for lag = 1:unit_lags, x_self(:,lag) = pred_spikes_shift(:,unit+(lag-1)*size(pred_spikes,2)); end
end


if ~isempty(kin_signals)
    if do_rcb
        kin_lags = unit_lags;
    end
    for lag = 1:length(kin_lags)
        if kin_lags(lag) == 0
            x_kin = [x_kin, cov_kin];
        else
            x_kin = [x_kin, cov_kin_shift(:,1+(kin_lags(lag)-1)*size(cov_kin,2):kin_lags(lag)*size(cov_kin,2))];
        end
    end
end

x_full = [x_self, x_kin, x_spikes];
x_basic = [x_self, x_kin];