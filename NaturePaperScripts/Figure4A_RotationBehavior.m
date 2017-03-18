%% Prep datasets and get paths etc
clear; clc; close all;
dataSummary;
trial_params;

% Session parameters
monkeys = {'Chewie','Mihili','all'};
tasks = {'CO'};
perts = {'VR'};
epochs = {'BL','AD','WO'};

session_idx = getFileDBidx(filedb, ...
    {'Task',tasks,'Perturbation',perts,'Monkey',monkeys}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});

idx_start = {'idx_go_cue', 0};
idx_end   = {'idx_trial_end', 0};

% Load data
fnames = cell(1,length(session_idx));
for iFile = 1:length(session_idx)
    file = session_idx(iFile);
    fnames{iFile} = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
end
func_calls = [trial_func_calls, {{@trimTD,idx_start,idx_end}}];
[trial_data,params] = loadTDfiles(fnames, func_calls{:});


%% MAKE ERROR PLOTS
% Get behavior metrics for all sessions
for monkey = monkeys
    if strcmpi(monkey{1},'all')
        idx = session_idx;
    else
        idx = session_idx(strcmpi(filedb.Monkey(session_idx),monkey{1}));
    end
    err = cell(length(idx),length(epochs));
    for iFile = 1:length(idx)
        file = idx(iFile);
        [~,td] = getTDidx(trial_data, ...
            'date',datestr(filedb.Date{file},'mm-dd-yyyy'), ...
            'monkey',filedb.Monkey{file}, ...
            'task',filedb.Task{file},...
            'perturbation',filedb.Perturbation{file});
        
        td = smoothSignals(td,struct('signals','vel','kernel_SD',0.1));
        temp = getLearningMetrics(td,struct( ...
            'time_window',{{'idx_go_cue',1;'idx_movement_on',16}}));
        for e = 1:length(epochs)
            if strcmpi(epochs{e},'AD')
                err_add = abs(td(1).perturbation_info(1));
            else
                err_add = 0;
            end
            flip_sign = sign(td(1).perturbation_info(1));
            err{iFile,e} = flip_sign*temp(getTDidx(td,'epoch',epochs{e})) + err_add;
        end
    end
    %err = cell2mat(reshape(cellfun(@(x) x(1:min(min(cellfun(@length,err),[],1)))',err,'uni',0),length(idx),1,length(epochs)));
    for e = 1:length(epochs)
        err(:,e) = cellfun(@(x) x(1:min(cellfun(@length,err(:,e)))),err(:,e),'uni',0);
    end
    
    figure;
    for e = 1:length(epochs)
        %m = squeeze(err(:,:,e))';
        m = cat(2,err{:,e});
        m = moving_average(m,2);
        
        subplot(1,length(epochs),e); hold all;
        %plot(squeeze(err(:,:,e))','-','Color',[0.4,0.4,0.4]);
        plot(m*180/pi,'-','Color',[0.6,0.6,0.6]);
        plot(mean(m,2)*180/pi,'k-','LineWidth',3);
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[-45 45]);
    end
end



%% MAKE PINWHEEL PLOTS
% get time-warped position traces for BL, early AD, late AD

monkey_files = [4,2];
num_reaches = 1; % first/last of this many to each target

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
    td = trimTD(td,{'idx_go_cue',5},{'idx_trial_end',-20});
    
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
        
        %         td_temp = trialAverage(td(getTDidx(td,'epoch','WO','target_direction',u(targ),'range',{'first',num_reaches})), ...
        %             'target_direction',struct('do_stretch',true,'num_samp',1000));
        %         plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'r-','LineWidth',2);
        %         td_temp = trialAverage(td(getTDidx(td,'epoch','WO','target_direction',u(targ),'range',{'last',num_reaches})), ...
        %             'target_direction',struct('do_stretch',true,'num_samp',1000));
        %         plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'r--','LineWidth',2);
    end
end



