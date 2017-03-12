%% Prep datasets and get paths etc
clear; clc; close all;
dataSummary;
trial_params;

% Session parameters
monkey = 'Chewie';
task = 'CO';
pert = 'FF';
date = '2016-10-07';

file = getFileDBidx(filedb, ...
    {'Task',task,'Perturbation',pert,'Monkey',monkey,'Date',date}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});

func_calls = [trial_func_calls, { ...
    {@removeBadNeurons,badneuron_params}, ...
    {@smoothSignals,struct('signals',{{'M1_spikes','PMd_spikes'}},'calc_fr',true,'kernel_SD',pn_kernel_SD)}, ...
    {@trimTD,idx_start,idx_end}}];

% load it
fname = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
[td,params] = loadTDfiles(fname, func_calls{:});

[~,td1] = getTDidx(td,'epoch','BL');
[~,td2] = getTDidx(td,'epoch','AD','range',[0.33 1]);
td_s = [td1, td2];

%% %%%%%%%% Plot stuff
% td = trimTD(td_s,{'idx_go_cue',-0},{'idx_go_cue',40});
td = appendTDs(trimTD(td_s,{'idx_go_cue',-0},{'idx_go_cue',40}),trimTD(td_s,{'idx_target_on',0},{'idx_target_on',50}));
td = binTD(td,5);
% td = softNormalize(td);
td = trialAverage(td,{'target_direction','epoch'});

% [~,td0] = getTDidx(td,'epoch','BL','range',[0 0.66]);
[~,td1] = getTDidx(td,'epoch','BL','range',[0 1]);
[~,td2] = getTDidx(td,'epoch','AD','range',[0.33 1]);

arrays = {'M1_spikes'};
params_corr = struct('signals',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
% [~,m1_bl_sort] = pairwiseCorr(td0,params_corr);
[m1_bl,m1_bl_sort] = pairwiseCorr(td1,params_corr);
[m1_ad,m1_ad_sort] = pairwiseCorr(td2,params_corr);
arrays = {'PMd_spikes'};
params_corr = struct('signals',{arrays},'cluster_order',true,'do_norm',true,'cluster_arrays',true);
% [~,pmd_bl_sort] = pairwiseCorr(td0,params_corr);
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

%% Now plot correlations against each other
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

temp = triu(m1_bl,1);
temp(temp == 0) = NaN;
m1_bl = temp;
temp = triu(pmd_bl,1);
temp(temp == 0) = NaN;
pmd_bl = temp;

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

