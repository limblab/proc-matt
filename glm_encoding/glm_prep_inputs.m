function varargout = glm_prep_inputs(trial_data,unit,trial_idx,varargin)
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
smooth_pca_proj_spikes = params.smooth_pca_proj_spikes;
center_pca = params.center_pca;
if ischar(do_pca)
    pca_dims = params.pca_dims.(cov_array); % get dimensions for the covariate array
end
do_rcb = params.do_rcb;
do_all_history = params.do_all_history;
do_self_history = params.do_self_history;
unit_lags = params.unit_lags;
kin_lags = params.kin_lags;
kin_signals = params.kin_signals;

% just trim it down to the trials we want
trial_data = trial_data(trial_idx);

%%% get the index of predicted unit
idx = true(1,size(trial_data(1).([cov_array '_spikes']),2));
if strcmpi(cov_array,pred_array)
    % take out this neuron and compute the new PCA
    idx(unit) = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get PCA. Gotta do it here if we want shift
if ischar(do_pca)
    for trial = 1:length(trial_data)
        temp = sqrt(trial_data(trial).([cov_array '_spikes']));
        if smooth_pca_proj_spikes
            temp = smoothSpikesForPCA(temp,params.bin_size,2*params.bin_size);
        end
        % de-mean data
        if center_pca
            temp = temp - repmat(params.pca_mu,size(temp,1),1);
        end
        temp = temp*params.pca_w;
        if ~ischar(pca_dims)
            temp = temp(:,pca_dims);
        end
        if strcmpi(do_pca,'potent')
            temp = temp * params.V_potent;
        elseif strcmpi(do_pca,'null')
            temp = temp * params.V_null;
        end
        trial_data(trial).([cov_array '_pca']) = temp;
    end
    if do_all_history
        trial_data = convBasisFunc(trial_data,[cov_array '_pca'],params);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get covariates
[cov_spikes, cov_spikes_shift, pred_spikes, pred_spikes_shift, cov_kin, cov_kin_shift] = get_covariates(trial_data,start_idx,end_idx,params);

[x_self, x_kin] = deal([]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% piece together covariates
if ischar(do_pca) % get PCA space
    %     % square root transform 4 lyfe
    %     x_spikes = sqrt(cov_spikes(:,idx))*params.pca_w;
    %     if ~ischar(params.pca_dims)
    %         x_spikes = x_spikes(:,pca_dims);
    %     end
    %
    %     % now project into potent/null spaces as desired
    %     if strcmpi(do_pca,'potent')
    %         x_spikes = x_spikes * params.V_potent;
    %     elseif strcmpi(do_pca,'null')
    %         x_spikes = x_spikes * params.V_null;
    %     end
    x_spikes = cov_spikes;
else % just use spike counts
    x_spikes = cov_spikes(:,idx);
end

% add shifted data as covariates if desired
if unit_lags > 0
    % add time shifted other neurons (NO PCA AT THE MOMENT)
    if do_all_history
        if do_pca
            x_spikes = [x_spikes, cov_spikes_shift];
        else
            idx_shift = repmat(idx,1,unit_lags);
            x_spikes = [x_spikes, cov_spikes_shift(:,idx_shift)];
        end
    end
    if do_self_history
        x_self = zeros(size(pred_spikes,1),unit_lags);
        for lag = 1:unit_lags, x_self(:,lag) = pred_spikes_shift(:,unit+(lag-1)*size(pred_spikes,2)); end
    end
end

if ~isempty(kin_signals)
    if do_rcb
        kin_lags = unit_lags;
        x_kin = [x_kin, cov_kin];
    end
    for lag = 1:length(kin_lags)
        if kin_lags(lag) == 0
            x_kin = [x_kin, cov_kin];
        else
            x_kin = [x_kin, cov_kin_shift(:,1+(kin_lags(lag)-1)*size(cov_kin,2):kin_lags(lag)*size(cov_kin,2))];
        end
    end
end

varargout{1} = pred_spikes(:,unit); % y
varargout{2} = [x_self, x_kin, x_spikes]; % x_full
varargout{3} = [x_self, x_kin]; % x_basic

