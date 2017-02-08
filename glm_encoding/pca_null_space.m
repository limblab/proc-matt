clear;
clc;

M1_ndim = 10;
PMd_ndim = 20;
num_samples = 2;

epochs = {'BL','AD','WO'};

% load trial_data
load('F:\TrialDataFiles\Chewie_CO_FF_2016-10-05.mat')


trial_data = truncateAndBin(trial_data,{'idx_target_on',0},{'idx_trial_end',0},num_samples);

%%
bl_idx = getTDidx(trial_data,'epoch','bl');
M1_spikes = cat(1,trial_data(bl_idx).M1_spikes);
PMd_spikes = cat(1,trial_data(bl_idx).PMd_spikes);

M1_spikes = M1_spikes(:,mean(M1_spikes,1) > 0);
PMd_spikes = PMd_spikes(:,mean(PMd_spikes,1) > 0);

M1_spikes = sqrt(M1_spikes);
PMd_spikes = sqrt(PMd_spikes);

% sort by highest firing rates and take first x PMd cells
% [~,I] = sort(sum(PMd_spikes,1));
% PMd_spikes = PMd_spikes(:,I(1:size(M1_spikes,2)));
temp = randperm(size(PMd_spikes,2));
PMd_spikes = PMd_spikes(:,temp(1:size(M1_spikes,2)));

%% now apply smoothing
disp('Smoothing the firing rates...');

bin_size            = 0.01*num_samples;
kernel_SD           = 3*bin_size;

binned_spikes_transf = M1_spikes;

% get nbr of channels and nbr of samples
[nbr_samples, nbr_chs]  = size(binned_spikes_transf);
% preallocate return matrix
smoothed_FR             = zeros(nbr_samples,nbr_chs);

% kernel half length is 3·SD out
kernel_hl               = ceil( 3 * kernel_SD / (bin_size) );
% create the kernel --it will have length 2*kernel_hl+1
kernel                  = normpdf( -kernel_hl*(bin_size) : ...
    bin_size : kernel_hl*(bin_size), ...
    0, kernel_SD );
% compute normalization factor --this factor depends on the number of taps
% actually used
nm                      = conv(kernel,ones(1,nbr_samples))';

% do the smoothing
for i = 1:nbr_chs
    aux_smoothed_FR     = conv(kernel,binned_spikes_transf(:,i)) ./ nm;
    % cut off the edges so that the result of conv is same length as the
    % original data
    smoothed_FR(:,i)    = aux_smoothed_FR(kernel_hl+1:end-kernel_hl);
end

M1_spikes = smoothed_FR;


binned_spikes_transf = PMd_spikes;

% get nbr of channels and nbr of samples
[nbr_samples, nbr_chs]  = size(binned_spikes_transf);
% preallocate return matrix
smoothed_FR             = zeros(nbr_samples,nbr_chs);

% kernel half length is 3·SD out
kernel_hl               = ceil( 3 * kernel_SD / (bin_size) );
% create the kernel --it will have length 2*kernel_hl+1
kernel                  = normpdf( -kernel_hl*(bin_size) : ...
    bin_size : kernel_hl*(bin_size), ...
    0, kernel_SD );
% compute normalization factor --this factor depends on the number of taps
% actually used
nm                      = conv(kernel,ones(1,nbr_samples))';

% do the smoothing
for i = 1:nbr_chs
    aux_smoothed_FR     = conv(kernel,binned_spikes_transf(:,i)) ./ nm;
    % cut off the edges so that the result of conv is same length as the
    % original data
    smoothed_FR(:,i)    = aux_smoothed_FR(kernel_hl+1:end-kernel_hl);
end

PMd_spikes = smoothed_FR;

%%
% build PCA model for M1
[M1_w, M1_scores, M1_eigen] = pca(M1_spikes);

% build PCA model for PMd
[PMd_w, PMd_scores, PMd_eigen] = pca(PMd_spikes);

%%
if 1
    figure,
    subplot(2,2,1);
    bar(M1_eigen/sum(M1_eigen));
    axis('tight');
    xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
    set(gca,'TickDir','out'),set(gca,'FontSize',14);
    xlim([0 size(M1_eigen,1)+1])
    
    subplot(2,2,3);
    plot(cumsum(M1_eigen/sum(M1_eigen)),'linewidth',3,'marker','d'),
    xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
    set(gca,'TickDir','out'),set(gca,'FontSize',14);
    xlim([0 size(M1_eigen,1)+1])
    ylim([0 1])
    
    subplot(2,2,2);
    bar(PMd_eigen/sum(PMd_eigen));
    axis('tight');
    xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
    set(gca,'TickDir','out'),set(gca,'FontSize',14);
    xlim([0 size(PMd_eigen,1)+1])
    
    subplot(2,2,4);
    plot(cumsum(PMd_eigen/sum(PMd_eigen)),'linewidth',3,'marker','d'),
    xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
    set(gca,'TickDir','out'),set(gca,'FontSize',14);
    xlim([0 size(PMd_eigen,1)+1])
    ylim([0 1])
end

%%
% build model from PMd to M1
y = M1_scores(:,1:M1_ndim);
x = [ones(size(PMd_scores,1),1), PMd_scores(:,1:PMd_ndim)];

% find the model
W = zeros( size(y,2), size(x,2) );
for i = 1:size(y,2)
    [b_pc, ~, ~, ~, stats_this] = regress(y(:,i),x);
    % fill MIMO matrix W
    W(i,:) = b_pc';
end


% do SVD of weights
[U, S, V]                   = svd( W );

% The output potent spaces is defined by the first m columns of V', where m
% is the number of dimensions of the output
V_potent                    = V(1:size(y,2),:)';
V_null                      = V(size(y,2)+1:end,:)';


%%
x_potent = x*V_potent;

