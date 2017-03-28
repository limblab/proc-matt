%% Prep datasets and get paths etc
clear; clc; close all;
dataSummary;

% Session parameters
monkey = 'Chewie';
task = 'CO';
pert = 'FF';
date = '2016-10-07';

file = getFileDBidx(filedb, ...
    {'Task',task,'Perturbation',pert,'Monkey',monkey,'Date',date}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});

idx_start = {'idx_target_on',0};
idx_end   = {'idx_trial_end', 0};

badtrial_params = struct(...
    'ranges', {{'idx_go_cue','idx_movement_on',[5 40]; ...
    'idx_peak_speed','idx_trial_end',[60 100]; ...
    'idx_target_on','idx_go_cue',[80 150]}});

badneuron_params_m1 = struct( ...
    'arrays','M1', ...
    'min_fr',1, ...
    'do_shunt_check',1, ...
    'use_trials',{{'epoch','BL'}});

badneuron_params_pmd = struct( ...
    'arrays','PMd', ...
    'min_fr',1, ...
    'do_shunt_check',1, ...
    'use_trials',{{'epoch','BL'}});

func_calls = { ...
    {@getTDidx,'result','R'}, ...
    {@removeBadTrials,badtrial_params}, ...
    {@removeBadNeurons,badneuron_params_m1}, ...
    {@removeBadNeurons,badneuron_params_pmd}, ...
    @sqrtTransform, ...
    {@smoothSignals,struct('signals',{{'M1_spikes','PMd_spikes'}},'calc_fr',true,'kernel_SD',0.1)}, ...
    {@trimTD,idx_start,idx_end}};

% load it
fname = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
[td,params] = loadTDfiles(fname, func_calls{:});

[~,td1] = getTDidx(td,'epoch','BL');
[~,td2] = getTDidx(td,'epoch','AD','range',[0.5 1]);
td_s = [td1, td2];

td = trimTD(td,{'idx_go_cue',-50},{'idx_go_cue',70});

% Find the potent and null spaces
pn_params = struct( ...
            'in_signals','PMd_spikes', ...
            'out_signals','M1_spikes', ...
            'in_dims',16, ...
            'out_dims',8, ...
            'use_trials',{{'epoch','BL'}}, ...
            'sqrt_transform',false, ...
            'do_smoothing',false);
td = getPotentSpace(td,pn_params);



td = trialAverage(td,{'epoch','target_direction'});
% td = subtractConditionMean(td);

%% Make plots
% first plot first two PMd PCs and first two M1 PCs
figure;
u = unique([td.target_direction]);
u = u([2,6]);

epochs = {'BL','AD'};
plot_colors = [0 0 1; 1 0 0];
plot_symbols = {'-','--'};

% find the max/min
ax_max = max(abs(cat(1,td.PMd_pca)));
ax_max = max(ax_max(1:3));
ax_min = min(cat(1,td.PMd_pca));
ax_min = min(ax_min(1:3));

plot_var = 'PMd_pca';
subplot(121); hold all;
for e = 1:length(epochs)
    for i = 1:length(u)
        td_idx = getTDidx(td,'epoch',epochs{e},'target_direction',u(i));
        for trial = td_idx
            temp = td(trial).(plot_var);
            plot3(temp(:,1),temp(:,2),temp(:,3),plot_symbols{e},'Color',plot_colors(i,:));
            plot3(temp(71,1),temp(71,2),temp(71,3),'o','Color',plot_colors(i,:));
        end
    end
end
axis('square');
set(gca,'Box','off','TickDir','out','XGrid','on','YGrid','on','ZGrid','on','FontSize',14);
% set(gca,'XLim',[ax_min, ax_max],'YLim',[ax_min, ax_max],'ZLim',[ax_min, ax_max]);
xlabel('PMd PC1');
ylabel('PMd PC2');
zlabel('PMd PC3');


subplot(122); hold all;
% find the max/min
ax_max = max(abs(cat(1,td.M1_pca)));
ax_max = max(ax_max(1:2));
ax_min = min(cat(1,td.M1_pca));
ax_min = min(ax_min(1:2));

plot_var = 'M1_pca';
subplot(122); hold all;
for e = 1:length(epochs)
    for i = 1:length(u)
        td_idx = getTDidx(td,'epoch',epochs{e},'target_direction',u(i));
        for trial = td_idx
            temp = td(trial).(plot_var);
            plot(temp(:,1),temp(:,2),plot_symbols{e},'Color',plot_colors(i,:));
            plot(temp(71,1),temp(71,2),'o','Color',plot_colors(i,:));
        end
    end
end
axis('square');
set(gca,'Box','off','TickDir','out','XGrid','on','YGrid','on','FontSize',14);
% set(gca,'XLim',[ax_min, ax_max],'YLim',[ax_min, ax_max]);
xlabel('M1 PC1');
ylabel('M1 PC2');


%%
epochs = {'BL','AD'};
plot_colors = [0 0 1; 1 0 0];
plot_symbols = {'-','--'};

u = unique([td.target_direction]);
u = u([2,6]);

figure;
subplot(121); hold all;
% plot first two potent and null
td_idx = getTDidx(td,'epoch','BL','target_direction',u(1));
for trial = td_idx
    plot(td(trial).PMdM1_potent(:,1),td(trial).PMdM1_null(:,1),'b-');
    plot(td(trial).PMdM1_potent(1,1),td(trial).PMdM1_null(1,1),'b+');
    plot(td(trial).PMdM1_potent(71,1),td(trial).PMdM1_null(71,1),'bo');
end
td_idx = getTDidx(td,'epoch','BL','target_direction',u(2));
for trial = td_idx
    plot(td(trial).PMdM1_potent(:,1),td(trial).PMdM1_null(:,1),'r-');
    plot(td(trial).PMdM1_potent(1,1),td(trial).PMdM1_null(1,1),'r+');
    plot(td(trial).PMdM1_potent(71,1),td(trial).PMdM1_null(71,1),'ro');
end
td_idx = getTDidx(td,'epoch','AD','target_direction',u(1));
for trial = td_idx
    plot(td(trial).PMdM1_potent(:,1),td(trial).PMdM1_null(:,1),'b--');
    plot(td(trial).PMdM1_potent(1,1),td(trial).PMdM1_null(1,1),'b+');
    plot(td(trial).PMdM1_potent(71,1),td(trial).PMdM1_null(71,1),'bo');
end
td_idx = getTDidx(td,'epoch','AD','target_direction',u(2));
for trial = td_idx
    plot(td(trial).PMdM1_potent(:,1),td(trial).PMdM1_null(:,1),'r--');
    plot(td(trial).PMdM1_potent(1,1),td(trial).PMdM1_null(1,1),'r+');
    plot(td(trial).PMdM1_potent(71,1),td(trial).PMdM1_null(71,1),'ro');
end
axis('square');
set(gca,'Box','off','TickDir','out','XGrid','on','YGrid','on','FontSize',14);
xlabel('Potent PC1');
ylabel('Null PC1');



subplot(122); hold all;
% find the max/min
ax_max = max(abs(cat(1,td.M1_pca)));
ax_max = max(ax_max(1:2));
ax_min = min(cat(1,td.M1_pca));
ax_min = min(ax_min(1:2));

plot_var = 'M1_pca';
for e = 1:length(epochs)
    for i = 1:length(u)
        td_idx = getTDidx(td,'epoch',epochs{e},'target_direction',u(i));
        for trial = td_idx
            temp = td(trial).(plot_var);
            plot(temp(:,1),temp(:,2),plot_symbols{e},'Color',plot_colors(i,:));
            plot(temp(1,1),temp(1,2),'+','Color',plot_colors(i,:));
            plot(temp(71,1),temp(71,2),'o','Color',plot_colors(i,:));
        end
    end
end
axis('square');
set(gca,'Box','off','TickDir','out','XGrid','on','YGrid','on','FontSize',14);
% set(gca,'XLim',[ax_min, ax_max],'YLim',[ax_min, ax_max]);
xlabel('M1 PC1');
ylabel('M1 PC2');

%% Plot M1/PMd against time
dim_idx = 1;

figure;
u = unique([td.target_direction]);
% u = u([2,6]);

epochs = {'BL','AD'};
plot_symbols = {'-','--'};

plot_colors = hex2rgb( ...
    {'327BA3',...
    '12AEB9',...
    'C0BC60',...
    'F2EA24',...
    'F8C337',...
    'E45C3B',...
    'A81E27',...
    '352F86'});

plot_var = 'PMd_pca';
subplot(411); hold all;
for e = 1:length(epochs)
    for i = 1:length(u)
        td_idx = getTDidx(td,'epoch',epochs{e},'target_direction',u(i));
        for trial = td_idx
            temp = td(trial).(plot_var);
            plot(temp(:,dim_idx),plot_symbols{e},'Color',plot_colors(i,:));
        end
    end
end
axis('tight')
V = axis;
plot([71 71],V(3:4),'k--');
set(gca,'Box','off','TickDir','out','FontSize',14);
xlabel('Time');
ylabel('PMd PC1');


plot_var = 'PMdM1_potent';
subplot(412); hold all;
for e = 1:length(epochs)
    for i = 1:length(u)
        td_idx = getTDidx(td,'epoch',epochs{e},'target_direction',u(i));
        for trial = td_idx
            temp = td(trial).(plot_var);
            plot(temp(:,dim_idx),plot_symbols{e},'Color',plot_colors(i,:));
        end
    end
end
axis('tight')
V = axis;
plot([71 71],V(3:4),'k--');
set(gca,'Box','off','TickDir','out','FontSize',14);
xlabel('Time');
ylabel('Null PC1');


plot_var = 'PMdM1_null';
subplot(413); hold all;
for e = 1:length(epochs)
    for i = 1:length(u)
        td_idx = getTDidx(td,'epoch',epochs{e},'target_direction',u(i));
        for trial = td_idx
            temp = td(trial).(plot_var);
            plot(temp(:,dim_idx),plot_symbols{e},'Color',plot_colors(i,:));
        end
    end
end
axis('tight')
V = axis;
plot([71 71],V(3:4),'k--');
set(gca,'Box','off','TickDir','out','FontSize',14);
xlabel('Time');
ylabel('Potent PC1');


subplot(414); hold all;
plot_var = 'M1_pca';
for e = 1:length(epochs)
    for i = 1:length(u)
        td_idx = getTDidx(td,'epoch',epochs{e},'target_direction',u(i));
        for trial = td_idx
            temp = td(trial).(plot_var);
            plot(temp(:,dim_idx),plot_symbols{e},'Color',plot_colors(i,:));
        end
    end
end
axis('tight')
V = axis;
plot([71 71],V(3:4),'k--');
set(gca,'Box','off','TickDir','out','FontSize',14);
xlabel('Time');
ylabel('M1 PC1');




