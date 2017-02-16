% load data
clear;
close all;
clc;

load('/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-07.mat');

[~,td1] = getTDidx(trial_data,'epoch','BL');
[~,td2] = getTDidx(trial_data,'epoch','AD','range',[0.66 1]);
td = [td1,td2];

td = truncateAndBin(td,1);
td = removeBadNeurons(td,struct('min_fr',5));
td = pruneBadTrials(td,struct('ranges', ...
    {{'idx_go_cue','idx_movement_on',[5,50]}}));

td = smoothSpikes(td,struct('sqrt_transform',true,'do_smoothing',true,'kernel_SD',0.1));
td = truncateAndBin(td,{'idx_movement_on',-50},{'idx_movement_on',40});
% td = truncateAndBin(td,{'idx_go_cue',0},{'idx_trial_end',-20});
td = trialAverage(td,{'target_direction','epoch'});

% td = softNormalize(td);


%%
% get subspaces
num_dims = 3;
pca_params = struct( ...
    'in_array','PMd', ...
    'out_array','M1', ...
    'in_dims',2*num_dims, ...
    'out_dims',num_dims, ...
    'do_smoothing',false, ...
    'sqrt_transform',false, ...
    'do_plot',false);

[td,temp] = getPotentSpace(td,pca_params);
td = getPCA(td,temp.w_out,temp.mu_out,struct('arrays','M1','do_smoothing',false,'sqrt_transform',false));

%%
close all;
clc;
figure; hold all;

dims = 1:10;%size(td(1).M1_spikes,2);
which_out = 'M1_pca';

which_one = 'PMd_pca';
[~,~,c_bl,b_bl] = linPredSignal(td(getTDidx(td,'epoch','BL')),which_one,which_out,num_dims);
[~,~,c_ad] = linPredSignal(td(getTDidx(td,'epoch','AD')),which_one,which_out,num_dims,b_bl);
subplot(131); hold all;
plot(c_bl(dims),'LineWidth',2)
plot(c_ad(dims),'LineWidth',2)
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[0,1]);

which_one = 'potent';
[~,~,c_bl,b_bl] = linPredSignal(td(getTDidx(td,'epoch','BL')),which_one,which_out,num_dims);
[~,~,c_ad] = linPredSignal(td(getTDidx(td,'epoch','AD')),which_one,which_out,num_dims,b_bl);
subplot(132); hold all;
plot(c_bl(dims),'LineWidth',2)
plot(c_ad(dims),'LineWidth',2)
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[0,1]);

which_one = 'null';
[~,~,c_bl,b_bl] = linPredSignal(td(getTDidx(td,'epoch','BL')),which_one,which_out,num_dims);
[~,~,c_ad] = linPredSignal(td(getTDidx(td,'epoch','AD')),which_one,which_out,num_dims,b_bl);
subplot(133); hold all;
plot(c_bl(dims),'LineWidth',2)
plot(c_ad(dims),'LineWidth',2)
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[0,1]);
