%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function varargout = getPCA(trial_data, varargin)
%
% [w, mu, eigen, scores, trial_data, params] = getPCA(trial_data, params);
%   Computes PCA projection for neural data. If you request trial_data as
% a final output, will add scores to each trial for use later. Must pass in
% 'array' field for struct. Note that this can be a cell of strings to pool
% data from multiple arrays. In this case, the rows of w will be as if you
% concatenated the two arrays together in the order they were provided.
%
% trial_data = getPCA(trial_data, w, mu, params);
%   Uses an old w and mu from previous getPCA call to add scores to
% trial_data as ARRAY_pca. Params is still needed to specify the array or
% if you want smoothing, etc.
%
%   NOTE: always de-means! Theoretically this could be modified to take a
% parameter to skip the de-meaning.
%
% INPUTS:
%   trial_data : the struct
%   params     : struct containing parameters
%     .array          : which units (can be cell array with multiple)
%     .trial_idx      : which trials to use (default: all)
%     .neurons        : which neurons to use (default: all) Note: for multiple arrays,
%                       neurons should be cell array with indices for each array
%     .sqrt_transform : flag to square root transform spikes (default: true)
%     .do_smoothing   : flag to convolve spikes with gaussian (default: true)
%     .bin_size       : size of time bins in trial_data (required for smoothing)
%     .kernel_SD      : kernel s.d. for smoothing (default: 2*bin_size)
%     .trial_avg      : flag to trial average (requires condition input) (default: false)
%     .trial_avg_cond : (string/cell) which conditions to average over
%     .do_plot        : flag to make scree plot (default: false)
%
% OUTPUTS:
%   w          : weight matrix for PCA projections
%   scores     : scores for the PCs
%   eigen      : eigenvalues for PC ranking
%   mu         : mean for each input (for de-meaning later)
%                  NOTE: if you use w later, you MUST demean using mu!!!
%   trial_data : old struct with added field for scores for each trial
%                   NOTE: if passing in old w and mu, only returns this
%
% EXAMPLES:
%   e.g. to compute covariance matrix
%       [w,mu,~,~,params] = getPCA(trial_data, struct('array','M1','bin_size',0.01));
%   e.g. to add scores to trial_data later using the above output
%       trial_data = getPCA(trial_data, w, mu, params);
% 
% Written by Matt Perich. Updated Feb 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = getPCA(trial_data, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
if length(varargin) == 1 % only params is provided
    new_pca = true;
    params = varargin{1};
else
    new_pca = false;
    if nargout > 1, error('When using old PCA, will only output trial_data'); end
    if length(varargin) == 3 % provided cov matrix and mu
        w = varargin{1};
        mu = varargin{2};
        params = varargin{3};
    else
        error('Incorrect number of inputs');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get defaults or process parameters
if ~exist('params','var'), params = struct(); end
if isfield(params,'array'),array = params.array; else, error('Need to specify array'); end
if isfield(params,'trial_idx'), trial_idx = params.trial_idx; else, trial_idx = 1:length(trial_data); end
if isfield(params,'neurons'), neurons = params.neurons; else, neurons = []; end
if isfield(params,'sqrt_transform'), sqrt_transform = params.sqrt_transform; else, sqrt_transform = true; end
if isfield(params,'do_smoothing'), do_smoothing = params.do_smoothing; else, do_smoothing = true; end
if isfield(params,'bin_size'), bin_size = params.bin_size; else, bin_size = NaN; end
if isfield(params,'kernel_SD'), kernel_SD = params.kernel_SD; else, kernel_SD = 2*bin_size; end
if isfield(params,'trial_avg'), trial_avg = params.trial_avg; else, trial_avg = false; end
if isfield(params,'trial_avg_cond'), trial_avg_cond = params.trial_avg_cond; else, trial_avg_cond = NaN; end
if isfield(params,'do_plot'), do_plot = params.do_plot; else, do_plot = false; end

if trial_avg && any(isnan(trial_avg_cond)), error('Must provide conditions to average trials over.'); end
if do_smoothing && isnan(bin_size), error('No bin size provided!'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process and prepare inputs
if ~iscell(array), array = {array}; end
if isempty(neurons)
    for i = 1:length(array), neurons{i} = 1:size(trial_data(1).([array{i} '_spikes']),2); end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop along trials to square root transform and smooth if desired
td = smoothSpikes(trial_data(trial_idx),params);
if trial_avg, td = trialAverage(td,trial_avg_cond); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concatenate specified trials
fr = [];
for i = 1:length(array)
    temp_fr = cat(1,td.([array{i} '_spikes']));
    fr = [fr, temp_fr(:,neurons{i})];
end
% get the time points that separate each trial later
trial_markers = [1,cumsum(cellfun(@(x) size(x,1),{td.pos}))];
trial_markers(end) = trial_markers(end)+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build PCA model for M1
if new_pca
    % compute PCA
    [w, scores, eigen,~,~,mu] = pca(fr,'Algorithm','svd');
    
    if do_plot
        figure,
        subplot(2,1,1);
        bar(eigen/sum(eigen));
        axis('tight');
        xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
        set(gca,'TickDir','out'),set(gca,'FontSize',14);
        xlim([0 size(fr,2)+1])
        
        subplot(2,1,2);
        plot(cumsum(eigen/sum(eigen)),'linewidth',3,'marker','d'),
        xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
        set(gca,'TickDir','out'),set(gca,'FontSize',14);
        xlim([0 size(fr,2)+1])
        ylim([0 1])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add scores to trial_data
if ~new_pca || nargout > 3
    if trial_avg, trial_data = td; end % THIS IS A HACK FOR NOW
    for trial = 1:length(trial_data)
        idx = trial_markers(trial):trial_markers(trial+1)-1;
        trial_data(trial).([[array{:}] '_pca']) = (fr(idx,:)-repmat(mu,length(idx),1))*w;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Package up outputs
if new_pca
    varargout{1} = w;
    varargout{2} = mu;
    varargout{3} = scores;
    varargout{4} = trial_data;
    varargout{5} = struct('array',array,'sqrt_transform',sqrt_transform,'do_smoothing',do_smoothing,'bin_size',bin_size,'kernel_SD',kernel_SD);
else
    varargout{1} = trial_data;
end