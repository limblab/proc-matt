%% Load data
clear;
clc;
close all;

filename = '/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-07.mat';
% Load data
load(filename);

trial_data = truncateAndBin(trial_data,{'idx_target_on',0},{'idx_trial_end',0});
[~,trial_data] = getTDidx(trial_data,'result','R');
[~,td1] = getTDidx(trial_data,'epoch','BL');
[~,td2] = getTDidx(trial_data,'epoch','AD','range',[0.66 1]);
td = [td1, td2];

td = smoothSignals(td,struct('signals',{getTDfields(td,'spikes')},'sqrt_transform',true,'calc_fr',true));
td_s = td;

%%
%%%%% PREPARATORY
td = truncateAndBin(td_s,{'idx_go_cue',-35},{'idx_go_cue',5});
td = trialAverage(td,{'target_direction','epoch'});

td = subtractConditionMean(td);

[~,td1] = getTDidx(td,'epoch','BL');
[~,td2] = getTDidx(td,'epoch','AD');

arrays = {'M1_spikes','PMd_spikes'};
params_corr = struct('signals',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
[rho_prep_bl,sort_prep_bl] = pairwiseCorr(td1,params_corr);
params_corr = struct('signals',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
[rho_prep_ad,sort_prep_ad] = pairwiseCorr(td2,params_corr);


%%%%%%%% MOVEMENT
td = truncateAndBin(td_s,{'idx_movement_on',5},{'idx_movement_on',35});
td = trialAverage(td,{'target_direction','epoch'});

[~,td1] = getTDidx(td,'epoch','BL');
[~,td2] = getTDidx(td,'epoch','AD');

params_corr = struct('signals',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
[rho_move_bl,sort_move_bl] = pairwiseCorr(td1,params_corr);
params_corr = struct('signals',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
[rho_move_ad,sort_move_ad] = pairwiseCorr(td2,params_corr);

%
sort1 = sort_prep_bl;

figure; subplot1(2,2);
subplot1(1); imagesc(flipud(rho_prep_bl(sort1,sort1))); axis('tight');
axis('square'); set(gca,'Box','off','TickDir','out','FontSize',14,'YTick',[],'XTick',[]);
title('Baseline','FontSize',14);
ylabel('Preparatory period','FontSize',14);
subplot1(2); imagesc(flipud(rho_prep_ad(sort1,sort1))); axis('tight');
axis('square'); set(gca,'Box','off','TickDir','out','FontSize',14,'YTick',[],'XTick',[]);
title('Learning','FontSize',14);

subplot1(3); imagesc(flipud(rho_move_bl(sort1,sort1))); axis('tight');
axis('square'); set(gca,'Box','off','TickDir','out','FontSize',14,'YTick',[],'XTick',[]);
ylabel('Movement period','FontSize',14);
subplot1(4); imagesc(flipud(rho_move_ad(sort1,sort1))); axis('tight');
axis('square'); set(gca,'Box','off','TickDir','out','FontSize',14,'YTick',[],'XTick',[]);

%%

%%%%%%%% MOVEMENT

td = truncateAndBin(td_s,{'idx_movement_on',-10},{'idx_movement_on',50});
td = removeBadNeurons(td,struct('min_fr',5));
td = softNormalize(td);
td = trialAverage(td,{'target_direction','epoch'});

[~,td0] = getTDidx(td,'epoch','BL','range',[0 0.66]);
[~,td1] = getTDidx(td,'epoch','BL','range',[0 1]);
[~,td2] = getTDidx(td,'epoch','AD','range',[0.5 1]);

arrays = {'M1_spikes'};
params_corr = struct('signals',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
[~,m1_bl_sort] = pairwiseCorr(td0,params_corr);
[m1_bl,m1_bl_sort] = pairwiseCorr(td1,params_corr);
[m1_ad,m1_ad_sort] = pairwiseCorr(td2,params_corr);
arrays = {'PMd_spikes'};
params_corr = struct('signals',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
[~,pmd_bl_sort] = pairwiseCorr(td0,params_corr);
[pmd_bl,pmd_bl_sort] = pairwiseCorr(td1,params_corr);
[pmd_ad,pmd_ad_sort] = pairwiseCorr(td2,params_corr);

figure; subplot1(2,2);
subplot1(1); imagesc(flipud(m1_bl(m1_bl_sort,m1_bl_sort))); axis('tight');
axis('square'); set(gca,'Box','off','TickDir','out','FontSize',14,'YTick',[],'XTick',[]);
title('M1','FontSize',14);
ylabel('Baseline','FontSize',14);

subplot1(2); imagesc(flipud(pmd_bl(pmd_bl_sort,pmd_bl_sort))); axis('tight');
axis('square'); set(gca,'Box','off','TickDir','out','FontSize',14,'YTick',[],'XTick',[]);
title('PMd','FontSize',14);

subplot1(3); imagesc(flipud(m1_ad(m1_bl_sort,m1_bl_sort))); axis('tight');
axis('square'); set(gca,'Box','off','TickDir','out','FontSize',14,'YTick',[],'XTick',[]);
ylabel('Adaptation','FontSize',14);

subplot1(4); imagesc(flipud(pmd_ad(pmd_bl_sort,pmd_bl_sort))); axis('tight');
axis('square'); set(gca,'Box','off','TickDir','out','FontSize',14,'YTick',[],'XTick',[]);

%% Plot pairwise in each condition
arrays = {'M1_spikes'};
params_corr = struct('signals',{arrays},'cluster_order',false,'do_norm',false,'cluster_arrays',false);
[m1_bl,m1_bl_sort] = pairwiseCorr(td1,params_corr);
[m1_ad,m1_ad_sort] = pairwiseCorr(td2,params_corr);
arrays = {'PMd_spikes'};
params_corr = struct('signals',{arrays},'cluster_order',false,'do_norm',false,'cluster_arrays',false);
[pmd_bl,pmd_bl_sort] = pairwiseCorr(td1,params_corr);
[pmd_ad,pmd_ad_sort] = pairwiseCorr(td2,params_corr);
arrays = {'M1_spikes','PMd_spikes'};
params_corr = struct('signals',{arrays},'cluster_order',false,'do_norm',false,'cluster_arrays',true);
[m1pmd_bl,m1pmd_bl_sort] = pairwiseCorr(td1,params_corr);
[m1pmd_ad,m1pmd_ad_sort] = pairwiseCorr(td2,params_corr);
m1pmd_bl = m1pmd_bl(1:size(m1_bl,1),size(m1_bl,1)+1:end);
m1pmd_ad = m1pmd_ad(1:size(m1_bl,1),size(m1_bl,1)+1:end);

figure;
subplot1(1,3);
subplot1(1); hold all;
plot([-1 1],[-1 1],'k--','LineWidth',2);
plot(reshape(m1_bl,1,numel(m1_bl)),reshape(m1_ad,1,numel(m1_ad)),'k.');
axis('square');
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-1 1],'YLim',[-1 1]);

subplot1(2); hold all;
plot([-1 1],[-1 1],'k--','LineWidth',2);
plot(reshape(pmd_bl,1,numel(pmd_bl)),reshape(pmd_ad,1,numel(pmd_ad)),'k.');
axis('square');
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-1 1],'YLim',[-1 1]);

subplot1(3); hold all;
plot([-1 1],[-1 1],'k--','LineWidth',2);
plot(reshape(m1pmd_bl,1,numel(m1pmd_bl)),reshape(m1pmd_ad,1,numel(m1pmd_ad)),'k.');
axis('square');
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-1 1],'YLim',[-1 1]);
