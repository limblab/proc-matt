clear;
clc;

M1_ndim = 10;
PMd_ndim = 20;
num_samples = 3;

epochs = {'BL','AD','WO'};

% load trial_data
load('F:\trial_data_files\Chewie_CO_FF_2016-09-19.mat')



%%
% adjust bin size
[M1_fr, PMd_fr] = deal([]);
epoch_indices = cell(1,length(epochs));
for idx_epoch = 1:length(epochs)
    [M1_fr_temp, PMd_fr_temp] = deal([]);
    td = trial_data(get_trial_data_indices(trial_data,'epoch',epochs{idx_epoch}));
    for idx_trial = 1:length(td)
        temp = full(trial_data(idx_trial).M1_spikes);
        bin_inds = 1:num_samples:size(temp,2);
        % fr is size bins x neurons
        fr = zeros(length(bin_inds),size(temp,1));
        for j = 1:length(bin_inds)-1
            fr(j,:) = sum(temp(:,bin_inds(j):bin_inds(j+1)),2);
        end
        M1_fr_temp = [M1_fr_temp; fr];
        
        temp = full(trial_data(idx_trial).PMd_spikes);
        % fr is size bins x neurons
        fr = zeros(length(bin_inds),size(temp,1));
        for j = 1:length(bin_inds)-1
            fr(j,:) = sum(temp(:,bin_inds(j):bin_inds(j+1)),2);
        end
        
        PMd_fr_temp = [PMd_fr_temp; fr];
        
    end, clear td temp bin_inds fr idx_trial;
    
    M1_fr = [M1_fr; M1_fr_temp];
    PMd_fr = [PMd_fr; PMd_fr_temp];
    epoch_indices{idx_epoch} = (1:size(M1_fr_temp,1))+sum(cellfun(@(x) length(x), epoch_indices));
end, clear M1_fr_temp PMd_fr_temp idx_epoch;

M1_fr = M1_fr(:,mean(M1_fr,1) > 0.1);
PMd_fr = PMd_fr(:,mean(PMd_fr,1) > 0.1);

M1_fr = sqrt(M1_fr);
PMd_fr = sqrt(PMd_fr);

%% now apply smoothing
disp('Smoothing the firing rates...');

    bin_size            = 0.01*num_samples;
    kernel_SD           = 2*bin_size;
    
    binned_spikes_transf = M1_fr;
    
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

M1_fr = smoothed_FR;


    binned_spikes_transf = PMd_fr;
    
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

PMd_fr = smoothed_FR;

%%
% build PCA model for M1
[M1_w, M1_scores, M1_eigen] = pca(M1_fr);

% build PCA model for PMd
[PMd_w, PMd_scores, PMd_eigen] = pca(PMd_fr);

%%
if 1
    figure,
    subplot(2,2,1);
    bar(M1_eigen/sum(M1_eigen));
    axis('tight');
    xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
    set(gca,'TickDir','out'),set(gca,'FontSize',14);
    xlim([0 size(M1_fr,2)+1])
    
    subplot(2,2,3);
    plot(cumsum(M1_eigen/sum(M1_eigen)),'linewidth',3,'marker','d'),
    xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
    set(gca,'TickDir','out'),set(gca,'FontSize',14);
    xlim([0 size(M1_fr,2)+1])
    ylim([0 1])
    
    subplot(2,2,2);
    bar(PMd_eigen/sum(PMd_eigen));
    axis('tight');
    xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
    set(gca,'TickDir','out'),set(gca,'FontSize',14);
    xlim([0 size(PMd_fr,2)+1])
    
    subplot(2,2,4);
    plot(cumsum(PMd_eigen/sum(PMd_eigen)),'linewidth',3,'marker','d'),
    xlabel('eigenvalue nbr.','FontSize',14),ylabel('explained variance','FontSize',14)
    set(gca,'TickDir','out'),set(gca,'FontSize',14);
    xlim([0 size(PMd_fr,2)+1])
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

%% build GLM model with PMd PCs
idx = epoch_indices{1};
idx_train = idx(1:floor(0.9*length(idx)));
idx_test = idx(floor(0.9*length(idx))+1:end);

%% use PMd spiking
x_spike = [ones(size(PMd_fr,1),1), PMd_fr];
pr2_spike = zeros(size(M1_fr,2),2);
b_spike = zeros(size(x_spike,2),size(M1_fr,2));
for i = 1:size(M1_fr,2)
    b_spike(:,i) = glmfit(x_spike(idx_train,2:end),M1_fr(idx_train,i),'poisson');
    
    yfit = exp([ones(length(idx_test),1), x_spike(idx_test,2:end)]*b_spike(:,i));
    
    try
        pr2_spike(i,:) = bootci(1000,{@compute_pseudo_R2,M1_fr(idx_test,i), yfit, mean(M1_fr(idx_test,i))},'Options',struct('UseParallel',true));
    catch
        pr2_spike(i,:) = NaN(2,1);
    end
    
    disp(i);    
end


%% use PMd PCs
x_pc = x;
pr2_pc = zeros(size(M1_fr,2),2);
b_pc = zeros(size(x_pc,2),size(M1_fr,2));
for i = 1:size(M1_fr,2)
    b_pc(:,i) = glmfit(x_pc(idx_train,2:end),M1_fr(idx_train,i),'poisson');
    
    yfit = exp([ones(length(idx_test),1), x_pc(idx_test,2:end)]*b_pc(:,i));
    
    pr2_pc(i,:) = bootci(1000,{@compute_pseudo_R2,M1_fr(idx_test,i), yfit, mean(M1_fr(idx_test,i))},'Options',struct('UseParallel',true));
    disp(i);
    
end

%% project PMd PCs into potent space and use as inputs to GLM
x_potent = x*V_potent;
pr2_potent = zeros(size(M1_fr,2),2);
b_potent = zeros(size(x_potent,2),size(M1_fr,2));
for i = 1:size(M1_fr,2)
    b_potent(:,i) = glmfit(x_potent(idx_train,2:end),M1_fr(idx_train,i),'poisson');
    
    yfit = exp([ones(length(idx_test),1), x_potent(idx_test,2:end)]*b_potent(:,i));
    
    pr2_potent(i,:) = bootci(1000,{@compute_pseudo_R2,M1_fr(idx_test,i), yfit, mean(M1_fr(idx_test,i))},'Options',struct('UseParallel',true));
    disp(i);
    
end


%% project PMd PCs into null space and use as inputs to GLM
x_null = x*V_null;
pr2_null = zeros(size(M1_fr,2),2);
b_null = zeros(size(x_null,2),size(M1_fr,2));
for i = 1:size(M1_fr,2)
    b_null(:,i) = glmfit(x_null(idx_train,2:end),M1_fr(idx_train,i),'poisson');
    
    yfit = exp([ones(length(idx_test),1), x_null(idx_test,2:end)]*b_null(:,i));
    
    pr2_null(i,:) = bootci(1000,{@compute_pseudo_R2,M1_fr(idx_test,i), yfit, mean(M1_fr(idx_test,i))},'Options',struct('UseParallel',true));
    disp(i);
    
end

%%
if 0
    figure;
    subplot(1,3,1); hist(pr2(:,1),-0.05:0.005:0.16); axis('tight'); set(gca,'Box','off','TickDir','out','FontSize',14); title('PCs','FontSize',16); ylabel('Count','FontSize',14); xlabel('pseudo-R2','FontSize',14);
    subplot(1,3,2); hist(pr2_potent(:,1),-0.05:0.005:0.16); axis('tight'); set(gca,'Box','off','TickDir','out','FontSize',14); title('Potent','FontSize',16);
    subplot(1,3,3); hist(pr2_null(:,1),-0.05:0.005:0.16); axis('tight'); set(gca,'Box','off','TickDir','out','FontSize',14); title('Null','FontSize',16);
end



%% now predict throughout AD

% cut out the crap
bins = 0:0.1:1;
idx = epoch_indices{2};
[pr2_ad_spike,pr2_ad_pc, pr2_ad_potent, pr2_ad_null] = deal(zeros(size(M1_fr,2),length(bins)-1,2));

for j = 2:length(bins)
    for i = 1:size(M1_fr,2)
        disp(['Bin ' num2str(j-1) ' neuron ' num2str(i)]);
        idx_test = idx(ceil(bins(j-1)*length(idx))+1:floor(bins(j)*length(idx)));
        
        
        yfit = exp([ones(length(idx_test),1), x_spike(idx_test,2:end)]*b_spike(:,i));
        try
            pr2_ad_spike(i,j-1,:) = bootci(1000,{@compute_pseudo_R2,M1_fr(idx_test,i), yfit, mean(M1_fr(idx_test,i))},'Options',struct('UseParallel',true));
        catch
            pr2_ad_spike(i,j-1,:) = NaN(1,2);
        end
        
        
        yfit = exp([ones(length(idx_test),1), x_pc(idx_test,2:end)]*b_pc(:,i));
        try
            pr2_ad_pc(i,j-1,:) = bootci(1000,{@compute_pseudo_R2,M1_fr(idx_test,i), yfit, mean(M1_fr(idx_test,i))},'Options',struct('UseParallel',true));
        catch
            pr2_ad_pc(i,j-1,:) = NaN(1,2);
        end
        
        
        yfit = exp([ones(length(idx_test),1), x_potent(idx_test,2:end)]*b_potent(:,i));
        try
            pr2_ad_potent(i,j-1,:) = bootci(1000,{@compute_pseudo_R2,M1_fr(idx_test,i), yfit, mean(M1_fr(idx_test,i))},'Options',struct('UseParallel',true));
        catch
            pr2_ad_potent(i,j-1,:) = NaN(1,2);
        end
        
        
        yfit = exp([ones(length(idx_test),1), x_null(idx_test,2:end)]*b_null(:,i));
        try
            pr2_ad_null(i,j-1,:) = bootci(1000,{@compute_pseudo_R2,M1_fr(idx_test,i), yfit, mean(M1_fr(idx_test,i))},'Options',struct('UseParallel',true));
        catch
            pr2_ad_null(i,j-1,:) = NaN(1,2);
        end
    end
end


%%
idx_spike =  pr2_spike(:,1) > 0 & ~any(isnan(pr2_ad_spike(:,:,1)),2);
idx_pc = pr2_pc(:,1) > 0 & ~any(isnan(pr2_ad_pc(:,:,1)),2);
idx_potent = pr2_potent(:,1) > 0 & ~any(isnan(pr2_ad_potent(:,:,1)),2);
idx_null = pr2_null(:,1) > 0 & ~any(isnan(pr2_ad_null(:,:,1)),2);

dpr2_spike = mean(pr2_ad_spike(idx_spike,:,:),3)-repmat(mean(pr2_spike(idx_spike,:),2),1,10);
dpr2_pc = mean(pr2_ad_pc(idx_pc,:,:),3)-repmat(mean(pr2_pc(idx_pc,:),2),1,10);
dpr2_potent = mean(pr2_ad_potent(idx_potent,:,:),3)-repmat(mean(pr2_potent(idx_potent,:),2),1,10);
dpr2_null = mean(pr2_ad_null(idx_null,:,:),3)-repmat(mean(pr2_null(idx_null,:),2),1,10);

dpr2_spike = dpr2_spike./repmat(mean(pr2_spike(idx_spike,:),2),1,10);
dpr2_pc = dpr2_pc./repmat(mean(pr2_pc(idx_pc,:),2),1,10);
dpr2_potent = dpr2_potent./repmat(mean(pr2_potent(idx_potent,:),2),1,10);
dpr2_null = dpr2_null./repmat(mean(pr2_null(idx_null,:),2),1,10);

m_spike = mean(dpr2_spike); s_spike = std(dpr2_spike)./sqrt(size(dpr2_spike,1));
m_pc = mean(dpr2_pc); s_pc = std(dpr2_pc)./sqrt(size(dpr2_pc,1));
m_potent = mean(dpr2_potent); s_potent = std(dpr2_potent)./sqrt(size(dpr2_potent,1));
m_null = mean(dpr2_null); s_null = std(dpr2_null)./sqrt(size(dpr2_null,1));


figure;
subplot(1,4,1); plot(dpr2_spike');
subplot(1,4,2); plot(dpr2_pc');
subplot(1,4,3); plot(dpr2_potent');
subplot(1,4,4); plot(dpr2_null');

figure;
subplot(1,4,1); hold all; plot(m_spike,'k'); plot(m_spike - s_spike,'k--'); plot(m_spike + s_spike,'k--');
subplot(1,4,2); plot(m_pc);
subplot(1,4,3); plot(m_potent);
subplot(1,4,4); plot(m_null);


