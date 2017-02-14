% load data
clear;
close all;
clc;

load('/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-07.mat');

%%
[~,td1] = getTDidx(trial_data,'epoch',{'BL'});
[~,td2] = getTDidx(trial_data,'epoch',{'AD'},'range',[0.5 1]);
trial_data = [td1 td2];

trial_data = truncateAndBin(trial_data,1);
trial_data = removeBadNeurons(trial_data,struct('min_fr',2));
trial_data = pruneBadTrials(trial_data,struct('ranges', ...
    {{'idx_go_cue','idx_movement_on',[10,40]}}));
%%
n_bins = 1;

align_idx = 'idx_target_on';
t_before = 20;
t_after = 30;
td1 = truncateAndBin(trial_data,n_bins,{align_idx,-t_before},{align_idx,t_after});

align_idx = 'idx_movement_on';
t_before = 20;
t_after = 50;
td2 = truncateAndBin(trial_data,n_bins,{align_idx,-t_before},{align_idx,t_after});

td = appendTDs(td1,td2);
%%
% get subspaces
pca_params = struct( ...
    'in_array','PMd', ...
    'out_array','M1', ...
    'in_dims',12, ...
    'out_dims',6, ...
    'sqrt_transform',true,...
    'do_smoothing', true, ...
    'bin_size', 0.01*n_bins, ...
    'kernel_SD',0.1, ...
    'trial_avg',true, ...
    'trial_avg_cond',{{'target_direction','epoch'}},...
    'do_plot',false);
[td,temp] = getPotentSpace(td,pca_params);
td = getPCA(td,temp.w_out,temp.mu_out,struct('array','M1','do_smoothing',false));

% figure;
% t_idx = 2:49;
% marker_loc = 20;
% u = unique([td.target_direction]);
% u = u(1:2:end);
% c = hsv(length(u));
% 
% which_one = 'null';
% 
% for i = 1:length(u)
%     idx =  getTDidx(td,'target_direction',u(i));
%     for j = 1:length(idx)
%         trial = idx(j);
%         hold all;
%         temp = td(trial).(which_one);
%         plot3(temp(t_idx,1),temp(t_idx,2),temp(t_idx,3),'-','LineWidth',3,'Color',c(i,:))
%         plot3(temp(t_idx(1),1),temp(t_idx(1),2),temp(t_idx(1),3),'o','LineWidth',3,'Color',c(i,:))
%         plot3(temp(marker_loc,1),temp(marker_loc,2),temp(marker_loc,3),'o','LineWidth',3,'Color',c(i,:))
%     end
% end
% title(which_one,'FontSize',14);
% axis('square');
% set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-0.6,0.6],'YLim',[-0.6,0.6],'ZLim',[-0.6,0.6]);
% xlabel('PC1','FontSize',14);
% ylabel('PC2','FontSize',14);
% zlabel('PC3','FontSize',14);
% grid on