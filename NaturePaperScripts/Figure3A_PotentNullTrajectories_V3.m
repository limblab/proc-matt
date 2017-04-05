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
    'min_fr',5, ...
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


% td = appendTDs(trimTD(td,{'idx_target_on',0},{'idx_target_on',50}),trimTD(td,{'idx_go_cue',0},{'idx_go_cue',70}));
td = trimTD(td,{'idx_go_cue',-60},{'idx_go_cue',80});

% Find the potent and null spaces
pn_params = struct( ...
    'in_signals','PMd_spikes', ...
    'out_signals','M1_spikes', ...
    'in_dims',4, ...
    'out_dims',2, ...
    'use_trials',{{'epoch',{'BL','AD'}}}, ...
    'sqrt_transform',false, ...
    'do_smoothing',false);
[td,pn_info] = getPotentSpace(td,pn_params);


% 
td = trialAverage(td,{'epoch','target_direction'});
% td = subtractConditionMean(td);


%% Plot M1/PMd against time
dim_idx = [1,2];

figure;
u = unique([td.target_direction]);
% u = u([2,6]);

epochs = {'BL'};
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
for d = 1:length(dim_idx)
    subplot(4,length(dim_idx),d); hold all;
    for e = 1:length(epochs)
        for i = 1:length(u)
            td_idx = getTDidx(td,'epoch',epochs{e},'target_direction',u(i));
            for trial = td_idx
                temp = td(trial).(plot_var);
                plot(temp(:,dim_idx(d)),plot_symbols{e},'Color',plot_colors(i,:));
            end
        end
    end
    axis('tight')
    V = axis;
    plot([61 61],V(3:4),'k--');
    set(gca,'Box','off','TickDir','out','FontSize',14);
    xlabel('Time');
    ylabel('PMd PC1');
    axis('square')
end


plot_var = 'PMdM1_null';
for d = 1:length(dim_idx)
    subplot(4,length(dim_idx),d+length(dim_idx)); hold all;
    for e = 1:length(epochs)
        for i = 1:length(u)
            td_idx = getTDidx(td,'epoch',epochs{e},'target_direction',u(i));
            for trial = td_idx
                temp = td(trial).(plot_var);
                plot(temp(:,dim_idx(d)),plot_symbols{e},'Color',plot_colors(i,:));
            end
        end
    end
    axis('tight')
    V = axis;
    plot([61 61],V(3:4),'k--');
    set(gca,'Box','off','TickDir','out','FontSize',14);
    xlabel('Time');
    ylabel('Null PC1');
    axis('square')
end

plot_var = 'PMdM1_potent';
for d = 1:length(dim_idx)
    subplot(4,length(dim_idx),2*length(dim_idx)+d); hold all;
    for e = 1:length(epochs)
        for i = 1:length(u)
            td_idx = getTDidx(td,'epoch',epochs{e},'target_direction',u(i));
            for trial = td_idx
                temp = td(trial).(plot_var);
                plot(temp(:,dim_idx(d)),plot_symbols{e},'Color',plot_colors(i,:));
            end
        end
    end
    axis('tight')
    V = axis;
    plot([61 61],V(3:4),'k--');
    set(gca,'Box','off','TickDir','out','FontSize',14);
    xlabel('Time');
    ylabel('Potent PC1');
    axis('square')
end

plot_var = 'M1_pca';
for d = 1:length(dim_idx)
    subplot(4,length(dim_idx),3*length(dim_idx)+d); hold all;
    for e = 1:length(epochs)
        for i = 1:length(u)
            td_idx = getTDidx(td,'epoch',epochs{e},'target_direction',u(i));
            for trial = td_idx
                temp = td(trial).(plot_var);
                plot(temp(:,dim_idx(d)),plot_symbols{e},'Color',plot_colors(i,:));
            end
        end
    end
    axis('tight')
    V = axis;
    plot([61 61],V(3:4),'k--');
    set(gca,'Box','off','TickDir','out','FontSize',14);
    xlabel('Time');
    ylabel('M1 PC1');
    axis('square')
end



