%% Prep datasets and get paths etc
clear; clc; close all;
dataSummary;
trial_params;

% Session parameters
monkeys = {'Chewie','Mihili','MrT'};
tasks = {'CO'};
perts = {'VR'};
epochs = {'BL','AD','WO'};

session_idx = getFileDBidx(filedb, ...
    {'Task',tasks,'Perturbation',perts,'Monkey',monkeys}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});

% Load data
fnames = cell(1,length(session_idx));
for iFile = 1:length(session_idx)
    file = session_idx(iFile);
    fnames{iFile} = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
end

func_calls = trial_func_calls;
[trial_data,params] = loadTDfiles(fnames, func_calls{:});


%% MAKE ERROR PLOTS
% Get behavior metrics for all sessions
for monkey = monkeys
    idx = session_idx(strcmpi(filedb.Monkey(session_idx),monkey{1}));
    err = cell(length(idx),length(epochs));
    for iFile = 1:length(idx)
        file = idx(iFile);
        [~,td] = getTDidx(trial_data, ...
            'date',datestr(filedb.Date{file},'mm-dd-yyyy'), ...
            'monkey',filedb.Monkey{file}, ...
            'task',filedb.Task{file},...
            'perturbation',filedb.Perturbation{file});
        td = smoothSignals(td,struct('signals','vel','kernel_SD',0.05));
        temp = getLearningMetrics(td,struct( ...
            'time_window',{{'idx_go_cue',0;'idx_movement_on',20}}));
        for e = 1:length(epochs)
            err{iFile,e} = temp(getTDidx(td,'epoch',epochs{e}));
        end
    end
    err = cell2mat(reshape(cellfun(@(x) x(1:min(min(cellfun(@length,err),[],1)))',err,'uni',0),length(idx),1,length(epochs)));
    figure;
    for e = 1:length(epochs)
        subplot(1,length(epochs),e);
        plot(squeeze(err(:,:,e))','.');
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14);
    end
end



%% MAKE PINWHEEL PLOTS
% get time-warped position traces for BL, early AD, late AD

monkey_files = [5,2,2];
num_reaches = 4; % first/last of this many to each target

for m = 1:length(monkeys)
    monkey = monkeys{m};
    
    figure; hold all;
    idx = session_idx(strcmpi(filedb.Monkey(session_idx),monkey));
    file = idx(monkey_files(m));
    
    [~,td] = getTDidx(trial_data, ...
        'date',datestr(filedb.Date{file},'mm-dd-yyyy'), ...
        'monkey',filedb.Monkey{file}, ...
        'task',filedb.Task{file},...
        'perturbation',filedb.Perturbation{file});
    td = trimTD(td,{'idx_go_cue',5},{'idx_trial_end',-10});
    
    % zero position
    pos_offsets = mean(cat(1,td.pos),1);
    
    u = unique([td.target_direction]);
    for targ = 1:length(u)
        td_temp = trialAverage(td(getTDidx(td,'epoch','BL','target_direction',u(targ))), ...
            'target_direction',struct('do_stretch',true,'num_samp',1000));
        plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'k-','LineWidth',2);
        
        td_temp = trialAverage(td(getTDidx(td,'epoch','AD','target_direction',u(targ),'range',{'first',num_reaches})), ...
            'target_direction',struct('do_stretch',true,'num_samp',1000));
        plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'b-','LineWidth',2);
        td_temp = trialAverage(td(getTDidx(td,'epoch','AD','target_direction',u(targ),'range',{'last',num_reaches})), ...
            'target_direction',struct('do_stretch',true,'num_samp',1000));
        plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'b--','LineWidth',2);
        
                td_temp = trialAverage(td(getTDidx(td,'epoch','WO','target_direction',u(targ),'range',{'first',num_reaches})), ...
            'target_direction',struct('do_stretch',true,'num_samp',1000));
        plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'r-','LineWidth',2);
        td_temp = trialAverage(td(getTDidx(td,'epoch','WO','target_direction',u(targ),'range',{'last',num_reaches})), ...
            'target_direction',struct('do_stretch',true,'num_samp',1000));
        plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'r--','LineWidth',2);
    end
end



