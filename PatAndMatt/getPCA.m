function [w, scores, eigen] = getPCA(trial_data, params)
% array can be cell array with multiple arrays. Then neurons should also be
% cell array with indices for each.
%
% Must provide array and bin size
if isfield(params,'array'),array = params.array; else error('Need to specify array'); end
if isfield(params,'bin_size'), bin_size = params.bin_size; else error('Must provide a bin size for smoothing'); end
if isfield(params,'trial_idx'), trial_idx = params.trial_idx; else trial_idx = 1:length(trial_data); end
if isfield(params,'neurons'), neurons = params.neurons; else neurons = []; end
if isfield(params,'kernel_SD'), kernel_SD = params.kernel_SD; else kernel_SD = 2*bin_size; end
if isfield(params,'do_plot'), do_plot = params.do_plot; else do_plot = false; end

% concatenate specified trials
if iscell(array) && length(array) > 1
    fr = [];
    for i = 1:length(array)
        if isempty(neurons), neurons{i} = 1:size(trial_data(1).([array{i} '_spikes']),2); end
        temp_fr = sqrt(cat(1,trial_data(trial_idx).([array{i} '_spikes'])));
        fr = [fr, temp_fr(:,neurons{i})];
    end
else
    if isempty(neurons), neurons = 1:size(trial_data(1).([array '_spikes']),2); end
    fr = sqrt(cat(1,trial_data(trial_idx).([array '_spikes'])));
    fr = fr(:,neurons);
end

% now apply smoothing
fr_smooth = smoothSpikesForPCA(fr,bin_size,kernel_SD);

% build PCA model for M1
[w, scores, eigen] = pca(fr_smooth);

if do_plot
    figure,
    subplot(2,1,1);
    bar(eigen/sum(eigen));
    axis('tight');
    xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
    set(gca,'TickDir','out'),set(gca,'FontSize',14);
    xlim([0 size(fr_smooth,2)+1])
    
    subplot(2,1,2);
    plot(cumsum(eigen/sum(eigen)),'linewidth',3,'marker','d'),
    xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
    set(gca,'TickDir','out'),set(gca,'FontSize',14);
    xlim([0 size(fr_smooth,2)+1])
    ylim([0 1])
end