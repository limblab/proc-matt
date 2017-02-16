%% Load data
clear;
clc;
close all;

filename = '/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-07.mat';
% Load data
load(filename);
trial_data = truncateAndBin(trial_data,{'idx_target_on',0},{'idx_trial_end',0});
[~,td1] = getTDidx(trial_data,'epoch','BL');
[~,td2] = getTDidx(trial_data,'epoch','AD','range',[0.5 1]);
td = [td1, td2];

td = pruneBadTrials(td,struct('ranges',{{'idx_go_cue','idx_movement_on',[5,50]}}));
td = removeBadNeurons(td,struct('min_fr',5));
td = smoothSpikes(td,struct('sqrt_transform',true,'do_smoothing',true,'calc_fr',true));
td = softNormalize(td);
td_s = td;

%%
%%%%% PREPARATORY
td = truncateAndBin(td_s,{'idx_go_cue',-35},{'idx_go_cue',5});
td = trialAverage(td,{'target_direction','epoch'});

td = subtractConditionMean(td);

[~,td1] = getTDidx(td,'epoch','BL');
[~,td2] = getTDidx(td,'epoch','AD');

arrays = {'M1','PMd'};
params_corr = struct('arrays',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
[rho_prep_bl,sort_prep_bl] = pairwiseCorr(td1,params_corr);
params_corr = struct('arrays',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
[rho_prep_ad,sort_prep_ad] = pairwiseCorr(td2,params_corr);


%%%%%%%% MOVEMENT
td = truncateAndBin(td_s,{'idx_movement_on',5},{'idx_movement_on',35});
td = trialAverage(td,{'target_direction','epoch'});

[~,td1] = getTDidx(td,'epoch','BL');
[~,td2] = getTDidx(td,'epoch','AD');

params_corr = struct('arrays',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
[rho_move_bl,sort_move_bl] = pairwiseCorr(td1,params_corr);
params_corr = struct('arrays',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
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
%%%%% PREPARATORY
td = truncateAndBin(td_s,{'idx_go_cue',-35},{'idx_go_cue',5});
td = trialAverage(td,{'target_direction','epoch'});

[~,td1] = getTDidx(td,'epoch','BL');
[~,td2] = getTDidx(td,'epoch','BL');

arrays = {'M1'};
params_corr = struct('arrays',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
[rho_prep_bl,sort_prep_bl] = pairwiseCorr(td1,params_corr);
arrays = {'PMd'};
params_corr = struct('arrays',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
[rho_prep_ad,sort_prep_ad] = pairwiseCorr(td2,params_corr);


%%%%%%%% MOVEMENT
td = truncateAndBin(td_s,{'idx_movement_on',5},{'idx_movement_on',35});
td = trialAverage(td,{'target_direction','epoch'});

[~,td1] = getTDidx(td,'epoch','BL');
[~,td2] = getTDidx(td,'epoch','AD');

arrays = {'M1'};
params_corr = struct('arrays',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
[rho_move_bl,sort_move_bl] = pairwiseCorr(td1,params_corr);
arrays = {'PMd'};
params_corr = struct('arrays',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
[rho_move_ad,sort_move_ad] = pairwiseCorr(td2,params_corr);

%
sort1 = sort_prep_bl;

figure; subplot1(2,2);
subplot1(1); imagesc(flipud(rho_prep_bl(sort1,sort1))); axis('tight');
axis('square'); set(gca,'Box','off','TickDir','out','FontSize',14,'YTick',[],'XTick',[]);
title('M1','FontSize',14);
ylabel('Preparatory period','FontSize',14);

sort1 = sort_prep_ad;
subplot1(2); imagesc(flipud(rho_prep_ad(sort1,sort1))); axis('tight');
axis('square'); set(gca,'Box','off','TickDir','out','FontSize',14,'YTick',[],'XTick',[]);
title('PMd','FontSize',14);

sort1 = sort_prep_bl;
subplot1(3); imagesc(flipud(rho_move_bl(sort1,sort1))); axis('tight');
axis('square'); set(gca,'Box','off','TickDir','out','FontSize',14,'YTick',[],'XTick',[]);
ylabel('Movement period','FontSize',14);
sort1 = sort_prep_ad;
subplot1(4); imagesc(flipud(rho_move_ad(sort1,sort1))); axis('tight');
axis('square'); set(gca,'Box','off','TickDir','out','FontSize',14,'YTick',[],'XTick',[]);