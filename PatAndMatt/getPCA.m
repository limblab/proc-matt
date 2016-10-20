function [w, scores, eigen] = getPCA(trial_data, params)
% Must provide array and bin size
if isfield(params,'array'),array = params.array; else error('Need to specify array'); end
if isfield(params,'bin_size'), bin_size = params.bin_size; else error('Must provide a bin size for smoothing'); end
if isfield(params,'trial_idx'), trial_idx = params.trial_idx; else trial_idx = 1:length(trial_data); end
if isfield(params,'neurons'), neurons = params.neurons; else neurons = 1:size(trial_data(1).([array '_spikes']),2); end
if isfield(params,'kernel_SD'), kernel_SD = params.kernel_SD; else kernel_SD = 2*bin_size; end

% concatenate specified trials
fr = sqrt(cat(1,trial_data(trial_idx).([array '_spikes'])));
fr = fr(:,neurons);

% now apply smoothing
% get nbr of channels and nbr of samples
[nbr_samples, nbr_chs]  = size(fr);
% preallocate return matrix
fr_smooth = zeros(nbr_samples,nbr_chs);

% kernel half length is 3·SD out
kernel_hl = ceil( 3 * kernel_SD / (bin_size) );
% create the kernel --it will have length 2*kernel_hl+1
kernel = normpdf( -kernel_hl*(bin_size) : ...
                            bin_size : kernel_hl*(bin_size), ...
                            0, kernel_SD );
% compute normalization factor --this factor depends on the number of taps
% actually used 
nm = conv(kernel,ones(1,nbr_samples))';

% do the smoothing
for i = 1:nbr_chs
    aux_smoothed_FR     = conv(kernel,fr(:,i)) ./ nm;
    % cut off the edges so that the result of conv is same length as the
    % original data
	fr_smooth(:,i)    = aux_smoothed_FR(kernel_hl+1:end-kernel_hl);
end

% build PCA model for M1
[w, scores, eigen] = pca(fr_smooth);

% for i = 1:length(trial_data)
%     trial_data(i).([array '_PCA']) = trial_data(i).([array '_spikes'])*w;
% end

%     figure,
%     subplot(2,2,1);
%     bar(M1_eigen/sum(M1_eigen));
%     axis('tight');
%     xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
%     set(gca,'TickDir','out'),set(gca,'FontSize',14);
%     xlim([0 size(M1_fr,2)+1])
%     
%     subplot(2,2,3);
%     plot(cumsum(M1_eigen/sum(M1_eigen)),'linewidth',3,'marker','d'),
%     xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
%     set(gca,'TickDir','out'),set(gca,'FontSize',14);
%     xlim([0 size(M1_fr,2)+1])
%     ylim([0 1])    