clear;
clc;
close all;

% plot population raster in BL and AD
load('F:\TrialDataFiles\Chewie_CO_FF_2016-09-15.mat')

% trial_data = truncateAndBin(trial_data,5);

%% Average across all trials to a target
close all;
figure; subplot1(3,2);

min_fr = 0.03;
targ_id = 1;

m1_means = mean(cat(1,trial_data.M1_spikes),1);
pmd_means = mean(cat(1,trial_data.PMd_spikes),1);

utheta = unique([trial_data.target_direction]);

bl_idx = find(getTDidx(trial_data,'epoch','bl','target_direction',utheta(targ_id)),1,'first');
bl_m1 = zeros(100,size(trial_data(1).M1_spikes,2),length(bl_idx));
bl_pmd = zeros(100,size(trial_data(1).PMd_spikes,2),length(bl_idx));
bl_vel = zeros(100,size(trial_data(1).vel,2),length(bl_idx));
for i = 1:length(bl_idx)
    trial = bl_idx(i);
    time_idx = trial_data(trial).idx_go_cue:trial_data(trial).idx_trial_end;
    
    temp = trial_data(trial).vel(time_idx,:);
    temp = smoothSpikesForPCA(temp,5,10);
    % interpolate to 100 points
    for j = 1:size(temp,2)
        bl_vel(:,j,i) = interp1(1:size(temp,1),temp(:,j),linspace(1,size(temp,1),100));
    end
    
    temp = trial_data(trial).M1_spikes(time_idx,:);
    temp = smoothSpikesForPCA(temp,5,5);
    temp = temp./repmat(m1_means,size(temp,1),1);
    for j = 1:size(temp,2)
        bl_m1(:,j,i) = interp1(1:size(temp,1),temp(:,j),linspace(1,size(temp,1),100));
    end
    
    temp = trial_data(trial).PMd_spikes(time_idx,:);
    temp = smoothSpikesForPCA(temp,5,5);
    temp = temp./repmat(pmd_means,size(temp,1),1);
    for j = 1:size(temp,2)
        bl_pmd(:,j,i) = interp1(1:size(temp,1),temp(:,j),linspace(1,size(temp,1),100));
    end
end

ad_idx = find(getTDidx(trial_data,'epoch','ad','target_direction',utheta(targ_id)),1,'last');
% ad_idx = ad_idx(floor(0.5*length(ad_idx)):end);

ad_m1 = zeros(100,size(trial_data(1).M1_spikes,2),length(ad_idx));
ad_pmd = zeros(100,size(trial_data(1).PMd_spikes,2),length(ad_idx));
ad_vel = zeros(100,size(trial_data(1).vel,2),length(ad_idx));
for i = 1:length(ad_idx)
    trial = ad_idx(i);
    time_idx = trial_data(trial).idx_go_cue:trial_data(trial).idx_trial_end;
    
    temp = trial_data(trial).vel(time_idx,:);
    temp = smoothSpikesForPCA(temp,5,10);
    % interpolate to 100 points
    for j = 1:size(temp,2)
        ad_vel(:,j,i) = interp1(1:size(temp,1),temp(:,j),linspace(1,size(temp,1),100));
    end
    
    temp = trial_data(trial).M1_spikes(time_idx,:);
    temp = smoothSpikesForPCA(temp,5,5);
    temp = temp./repmat(m1_means,size(temp,1),1);
    for j = 1:size(temp,2)
        ad_m1(:,j,i) = interp1(1:size(temp,1),temp(:,j),linspace(1,size(temp,1),100));
    end
    
    temp = trial_data(trial).PMd_spikes(time_idx,:);
    temp = smoothSpikesForPCA(temp,5,5);
    temp = temp./repmat(pmd_means,size(temp,1),1);
    for j = 1:size(temp,2)
        ad_pmd(:,j,i) = interp1(1:size(temp,1),temp(:,j),linspace(1,size(temp,1),100));
    end
end

max_vel = max([max(max(mean(bl_vel,3))), max(max(mean(ad_vel,3)))]);
min_vel = min([min(min(mean(bl_vel,3))), min(min(mean(ad_vel,3)))]);

max_m1 = max([max(max(mean(bl_m1,3))), max(max(mean(ad_m1,3)))]);
min_m1 = min([min(min(mean(bl_m1,3))), min(min(mean(ad_m1,3)))]);

max_pmd = max([max(max(mean(bl_pmd,3))), max(max(mean(ad_pmd,3)))]);
min_pmd = min([min(min(mean(bl_pmd,3))), min(min(mean(ad_pmd,3)))]);

subplot1(1);
plot(mean(bl_vel,3),'LineWidth',2);
axis('tight');
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[min_vel, max_vel]);
subplot1(3);
good_idx = mean(cat(1,trial_data.M1_spikes),1) > min_fr;
pds = getNeuronTuning(trial_data(getTDidx(trial_data,'epoch','bl')),'regress',struct('window',{{'idx_movement_on',-30; 'idx_trial_end',0}},'array','M1'));
[~,I] = sort(pds);
imagesc(mean(bl_m1(:,good_idx(I),:),3)');
axis('tight');
colormap('gray');
set(gca,'Box','off','TickDir','out','FontSize',14);
subplot1(5);
good_idx = mean(cat(1,trial_data.PMd_spikes),1) > min_fr;
pds = getNeuronTuning(trial_data(getTDidx(trial_data,'epoch','bl')),'regress',struct('window',{{'idx_movement_on',-30; 'idx_trial_end',0}},'array','PMd'));
[~,I] = sort(pds);
imagesc(mean(bl_pmd(:,good_idx(I),:),3)');
axis('tight');
colormap('gray');
set(gca,'Box','off','TickDir','out','FontSize',14);

subplot1(2);
plot(mean(ad_vel,3),'LineWidth',2);
axis('tight');
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[min_vel, max_vel]);
subplot1(4);
good_idx = mean(cat(1,trial_data.M1_spikes),1) > min_fr;
pds = getNeuronTuning(trial_data(getTDidx(trial_data,'epoch','bl')),'regress',struct('window',{{'idx_movement_on',-30; 'idx_trial_end',0}},'array','M1'));
[~,I] = sort(pds);
imagesc(mean(ad_m1(:,good_idx(I),:),3)');
axis('tight');
colormap('gray');
set(gca,'Box','off','TickDir','out','FontSize',14);
subplot1(6);
good_idx = mean(cat(1,trial_data.PMd_spikes),1) > min_fr;
pds = getNeuronTuning(trial_data(getTDidx(trial_data,'epoch','bl')),'regress',struct('window',{{'idx_movement_on',-30; 'idx_trial_end',0}},'array','PMd'));
[~,I] = sort(pds);
imagesc(mean(ad_pmd(:,good_idx(I),:),3)');
axis('tight');
colormap('gray');
set(gca,'Box','off','TickDir','out','FontSize',14);


