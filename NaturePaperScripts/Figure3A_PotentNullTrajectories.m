function Figure3A_PotentNullTrajectories()
clear; clc; close all;
dataSummary;
trial_params;

% Session parameters
monkey = 'Chewie';
task = 'CO';
pert = 'FF';
date = '2016-10-07';

pn_params = struct( ...
    'in_signals','PMd_spikes', ...
    'out_signals','M1_spikes', ...
    'in_dims',10, ...
    'out_dims',5, ...
    'use_trials',{{'epoch','BL'}}, ...
    'trim_idx',{{'idx_go_cue',-50;'idx_go_cue',50}});

file = getFileDBidx(filedb, ...
    {'Task',task,'Perturbation',pert,'Monkey',monkey,'Date',date}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});

func_calls = [trial_func_calls, { ...
    {@removeBadNeurons,badneuron_params}, ...
    {@getTDidx,'epoch',{'BL','AD'}}, ...
    {@smoothSignals,struct('signals',{{'M1_spikes','PMd_spikes'}},'calc_fr',true,'kernel_SD',pn_kernel_SD)}, ...
    {@trimTD,idx_start,idx_end}, ...
    {@getPotentSpace,pn_params}}];

% load it
fname = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
[trial_data,params] = loadTDfiles(fname, func_calls{:});

[~,td1] = getTDidx(trial_data,'epoch',{'BL'});
[~,td2] = getTDidx(trial_data,'epoch',{'AD'},'range',[0.5 1]);
trial_data = [td1 td2];


align_idx = 'idx_target_on';
t_before = 10;
t_after = 30;
td1 = trimTD(trial_data,{align_idx,-t_before},{align_idx,t_after});
align_idx = 'idx_movement_on';
t_before = 10;
t_after = 50;
td2 = trimTD(trial_data,{align_idx,-t_before},{align_idx,t_after});
td = appendTDs(td1,td2);

td = trialAverage(td,{'target_direction','epoch'});

potentspace_plot_script(td);
end

%%%%%%%%%%%%%
function potentspace_plot_script(td)
close all;
figure;

epochs = {'BL','AD'};
arrow_size = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_idx = 2:50;
marker_loc = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputArgs = struct( ...
    'which_one', 'PMdM1_null', ...
    'axis_range',[-1, 1], ...
    'sp_num',1, ...
    't_idx',t_idx, ...
    'marker_loc',marker_loc, ...
    'arrow_size',arrow_size);
do_the_plot(td,epochs,inputArgs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputArgs = struct( ...
    'which_one', 'PMdM1_potent', ...
    'axis_range',[-1, 1], ...
    'sp_num',2, ...
    't_idx',t_idx, ...
    'marker_loc',marker_loc, ...
    'arrow_size',arrow_size);
do_the_plot(td,epochs,inputArgs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputArgs = struct( ...
    'which_one', 'M1_pca', ...
    'axis_range',[-1, 1.4], ...
    'sp_num',3, ...
    't_idx',t_idx, ...
    'marker_loc',marker_loc, ...
    'arrow_size',arrow_size+2);
do_the_plot(td,epochs,inputArgs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_idx = 55:80;
marker_loc = 70;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputArgs = struct( ...
    'which_one', 'PMdM1_null', ...
    'axis_range',[-1, 1], ...
    'sp_num',1, ...
    't_idx',t_idx, ...
    'marker_loc',marker_loc, ...
    'arrow_size',arrow_size);
do_the_plot(td,epochs,inputArgs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputArgs = struct( ...
    'which_one', 'PMdM1_potent', ...
    'axis_range',[-1, 1], ...
    'sp_num',2, ...
    't_idx',t_idx, ...
    'marker_loc',marker_loc, ...
    'arrow_size',4);
do_the_plot(td,epochs,inputArgs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputArgs = struct( ...
    'which_one', 'M1_pca', ...
    'axis_range',[-1, 1.4], ...
    'sp_num',3, ...
    't_idx',t_idx, ...
    'marker_loc',marker_loc, ...
    'arrow_size',arrow_size+2);
do_the_plot(td,epochs,inputArgs);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function do_the_plot(td,epochs,inputArgs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
which_one = '';
axis_range = [];
sp_num = 0;
t_idx = [];
marker_loc = 0;
arrow_size = 0;
assignParams(who,inputArgs); %get parameters

sp_cols = 3;
sp_rows = 1;

epoch_sym = {'-','--'};
epoch_widths = [1,2];
axis_min = axis_range(1);
axis_max = axis_range(2);

u = unique([td.target_direction]);
c = phasemap(length(u));

subplot(sp_rows,sp_cols,sp_num);
for i = 1:length(u)
    for e = 1:length(epochs)
        idx =  getTDidx(td,'target_direction',u(i),'epoch',epochs{e});
        for j = 1:length(idx)
            trial = idx(j);
            hold all;
            temp = td(trial).(which_one);
            plot3(temp(t_idx,1),temp(t_idx,2),temp(t_idx,3),epoch_sym{e},'LineWidth',epoch_widths(e),'Color',c(i,:))
            plot3(temp(marker_loc,1),temp(marker_loc,2),temp(marker_loc,3),'o','LineWidth',1,'Color',c(i,:))
            arrowMMC3(temp(t_idx(end-2),1:3),temp(t_idx(end-1),1:3),temp(t_idx(end),1:3),arrow_size, ...
                [axis_min,axis_max,axis_min,axis_max,axis_min,axis_max],c(i,:),c(i,:));
        end
    end
end
title(which_one,'FontSize',14);
axis('square');

set(gca,'Box','off','TickDir','out','FontSize',14);
xlabel('PC1','FontSize',14);
ylabel('PC2','FontSize',14);
zlabel('PC3','FontSize',14);
grid on

end