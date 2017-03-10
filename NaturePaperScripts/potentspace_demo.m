% load data
clear;
close all;
clc;

load('/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-05.mat');

%%
[~,trial_data] = getTDidx(trial_data,'result','R');
[~,td1] = getTDidx(trial_data,'epoch',{'BL'});
[~,td2] = getTDidx(trial_data,'epoch',{'AD'},'range',[0.5 1]);
trial_data = [td1 td2];

trial_data = removeBadNeurons(trial_data,struct('min_fr',2));
trial_data = removeBadTrials(trial_data,struct('ranges', ...
    {{'idx_go_cue','idx_movement_on',[10,40]}}));
%%

align_idx = 'idx_target_on';
t_before = 10;
t_after = 30;
td1 = trimTD(trial_data,{align_idx,-t_before},{align_idx,t_after});

align_idx = 'idx_movement_on';
t_before = 10;
t_after = 50;
td2 = trimTD(trial_data,{align_idx,-t_before},{align_idx,t_after});

td = appendTDs(td1,td2);

td = sqrtTransform(td);
td = smoothSignals(td,struct('signals',{{'M1_spikes','PMd_spikes'}},'kernel_SD',0.1));
td = trialAverage(td,{'target_direction','epoch'});

% td = subtractConditionMean(td);

% get subspaces
pca_params = struct( ...
    'in_signals','PMd_spikes', ...
    'out_signals','M1_spikes', ...
    'in_dims',12, ...
    'out_dims',6, ...
    'sqrt_transform',false,...
    'do_smoothing', false, ...
    'kernel_SD',0.1, ...
    'do_plot',false);
[td,temp] = getPotentSpace(td,pca_params);



%%

align_idx = 'idx_movement_on';
t_before = 50;
t_after = 50;
td = trimTD(trial_data,{align_idx,-t_before},{align_idx,t_after});

td = sqrtTransform(td);
td = smoothSignals(td,struct('signals',{{'M1_spikes','PMd_spikes'}},'kernel_SD',0.1));

% td = subtractConditionMean(td);

% get subspaces
pca_params = struct( ...
    'in_signals','PMd_spikes', ...
    'out_signals','M1_spikes', ...
    'in_dims',12, ...
    'out_dims',8, ...
    'kernel_SD',0.1, ...
    'do_plot',false);
[td,temp] = getPotentSpace(td,pca_params);

td = trialAverage(td,{'target_direction','epoch'});

figure;
t_idx = 1:100;
marker_loc = 51;
u = unique([td.target_direction]);
c = phasemap(length(u));

for i = 1:length(u)
    idx =  getTDidx(td,'target_direction',u(i),'epoch','BL');
    for j = 1:length(idx)
        trial = idx(j);
        hold all;
        
        temp1 = td(trial).PMdM1_null;
        temp2 = td(trial).PMdM1_potent;
        temp = [temp2(:,1:2),temp1(:,1)];
        
        plot3(temp(t_idx,1),temp(t_idx,2),temp(t_idx,3),'-','LineWidth',3,'Color',c(i,:))
        plot3(temp(t_idx(1),1),temp(t_idx(1),2),temp(t_idx(1),3),'o','LineWidth',3,'Color',c(i,:))
        plot3(temp(marker_loc,1),temp(marker_loc,2),temp(marker_loc,3),'o','LineWidth',3,'Color',c(i,:))
    end
end
title('BL','FontSize',14);
axis('square');
set(gca,'Box','off','TickDir','out','FontSize',14);
xlabel('PC1','FontSize',14);
ylabel('PC2','FontSize',14);
zlabel('PC3','FontSize',14);
grid on
