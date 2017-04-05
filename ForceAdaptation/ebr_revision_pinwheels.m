%% Prep datasets and get paths etc
clear; clc; close all;
dataSummary;
trial_params;

% Session parameters
monkeys = {'Chewie','Mihili'};
tasks = {'CO'};
perts = {'FF'};
epochs = {'BL','AD','WO'};

session_idx = getFileDBidx(filedb, ...
    {'Task',tasks,'Perturbation',perts,'Monkey',monkeys}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-06-16''))', ...
    '~(ismember(filedb.Monkey,''Chewie'') & datenum(filedb.Date) > datenum(''2015-07-08''))'});

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

%% MAKE PINWHEEL PLOTS
% get time-warped position traces for BL, early AD, late AD
close all;
monkeys = {'Mihili'};
monkey_files = [2];
num_reaches = 1; % first/last of this many to each target

for m = 1:length(monkeys)
    monkey = monkeys{m};
    
    figure; hold all;
    idx = session_idx(strcmpi(filedb.Monkey(session_idx),monkey));
    %file = idx(monkey_files(m));
    file = idx(monkey_files(m));
    
    [~,td] = getTDidx(trial_data, ...
        'date',datestr(filedb.Date{file},'mm-dd-yyyy'), ...
        'monkey',filedb.Monkey{file}, ...
        'task',filedb.Task{file},...
        'perturbation',filedb.Perturbation{file});
    td = trimTD(td,{'idx_go_cue',0},{'idx_trial_end',-0});
    
    % zero position
    pos_offsets = mean(cat(1,td.pos),1);
    if td(1).perturbation_info(2) < 0
        for trial = 1:length(td)
            td(trial).pos = -td(trial).pos;
        end
    end
    u = unique([td.target_direction]);
    for targ = 1:length(u)
        td_reach = td(getTDidx(td,'epoch','BL','target_direction',u(targ),'range',{'last',num_reaches}));
        for j = 1:num_reaches
            td_temp = trialAverage(td_reach(j),'target_direction',struct('do_stretch',true,'num_samp',1000));
            plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'k-','LineWidth',2);
        end
        
        td_reach = td(getTDidx(td,'epoch','AD','target_direction',u(targ),'range',{'first',num_reaches}));
        for j = 1:num_reaches
            td_temp = trialAverage(td_reach(j),'target_direction',struct('do_stretch',true,'num_samp',1000));
            plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'b-','LineWidth',2);
        end
        
        td_reach = td(getTDidx(td,'epoch','AD','target_direction',u(targ),'range',{'last',num_reaches}));
        for j = 1:num_reaches
            td_temp = trialAverage(td_reach(j),'target_direction',struct('do_stretch',true,'num_samp',1000));
            plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'b--','LineWidth',2);
        end
        
        td_reach = td(getTDidx(td,'epoch','WO','target_direction',u(targ),'range',{'first',num_reaches}));
        for j = 1:num_reaches
            td_temp = trialAverage(td_reach(j),'target_direction',struct('do_stretch',true,'num_samp',1000));
            plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'r-','LineWidth',2);
        end
    end
    legend({'BL','First AD','Last AD','First WO'});
    axis('square');
    set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-9 9],'YLim',[-9,9]);
    
end


%% MAKE PINWHEEL PLOTS
% get time-warped position traces for BL, early AD, late AD
close all;
monkeys = {'Chewie','Mihili'};
num_reaches = 1; % first/last of this many to each target

for m = 1:length(monkeys)
    monkey = monkeys{m};
    
    figure; hold all;
    idx = session_idx(strcmpi(filedb.Monkey(session_idx),monkey));
    for file = idx'
        [~,td] = getTDidx(trial_data, ...
            'date',datestr(filedb.Date{file},'mm-dd-yyyy'), ...
            'monkey',filedb.Monkey{file}, ...
            'task',filedb.Task{file},...
            'perturbation',filedb.Perturbation{file});
        td = trimTD(td,{'idx_go_cue',0},{'idx_trial_end',-0});
        
        if td(1).perturbation_info(2) < 0
        % zero position
        pos_offsets = mean(cat(1,td.pos),1);
%         if td(1).perturbation_info(2) < 0
%             for trial = 1:length(td)
%                 td(trial).pos = -td(trial).pos;
%             end
%         end
        u = unique([td.target_direction]);
        for targ = 1:length(u)
            td_temp = trialAverage(td(getTDidx(td,'epoch','BL','target_direction',u(targ),'range',{'last',num_reaches})), ...
                'target_direction',struct('do_stretch',true,'num_samp',100));
            plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'k-');
            
            td_temp = trialAverage(td(getTDidx(td,'epoch','AD','target_direction',u(targ),'range',{'first',num_reaches})), ...
                'target_direction',struct('do_stretch',true,'num_samp',100));
            plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'b-');
            td_temp = trialAverage(td(getTDidx(td,'epoch','AD','target_direction',u(targ),'range',{'last',num_reaches})), ...
                'target_direction',struct('do_stretch',true,'num_samp',100));
            plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'b--');
            
            td_temp = trialAverage(td(getTDidx(td,'epoch','WO','target_direction',u(targ),'range',{'first',num_reaches})), ...
                'target_direction',struct('do_stretch',true,'num_samp',100));
            plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'r-');
            %                 td_temp = trialAverage(td(getTDidx(td,'epoch','WO','target_direction',u(targ),'range',{'last',num_reaches})), ...
            %                     'target_direction',struct('do_stretch',true,'num_samp',1000));
            %                 plot(td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2),'r--','LineWidth',2);
        end
        end
%         legend({'BL','First AD','Last AD','First WO'});
    end
    axis('square');
    set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-10 10],'YLim',[-10,10]);
    
end


%% just do first reach and rotate
% get time-warped position traces for BL, early AD, late AD
close all;
monkeys = {'Chewie','Mihili'};
num_reaches = 1; % first/last of this many to each target

for m = 1:length(monkeys)
    monkey = monkeys{m};
    
    figure; hold all;
    idx = session_idx(strcmpi(filedb.Monkey(session_idx),monkey));
    %file = idx(monkey_files(m));
    for file = idx'
        [~,td] = getTDidx(trial_data, ...
            'date',datestr(filedb.Date{file},'mm-dd-yyyy'), ...
            'monkey',filedb.Monkey{file}, ...
            'task',filedb.Task{file},...
            'perturbation',filedb.Perturbation{file});
        td = trimTD(td,{'idx_go_cue',5},{'idx_trial_end',-20});
        
        % zero position
        pos_offsets = mean(cat(1,td.pos),1);
        
        
                %%%%%%%%%%%%%%
        td_temp = trialAverage(td(getTDidx(td,'epoch','AD','range',{'first',num_reaches})), ...
            'target_direction',struct('do_stretch',true,'num_samp',1000));
        pos = [td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2)];
        ad_dir = td_temp.target_direction;
        % rotate position by
        R = [cos(-td_temp.target_direction) -sin(-td_temp.target_direction); sin(-td_temp.target_direction) cos(-td_temp.target_direction)];
        newPos = zeros(size(pos));
        for j = 1:length(pos)
            newPos(j,:) = R*(pos(j,:)');
        end
        newPos = zeros(size(pos));
        for j = 1:length(pos)
            newPos(j,:) = R*(pos(j,:)');
        end
        if td(1).perturbation_info(2) < 0
            newPos = -newPos;
        end
        newPos = newPos - repmat(newPos(1,:),size(newPos,1),1);
        plot(newPos(:,1),newPos(:,2),'b-','LineWidth',2);
        
        
        
        %%%%%%%
        td_temp = trialAverage(td(getTDidx(td,'epoch','BL','target_direction',ad_dir,'range',{'last',num_reaches})), ...
            'target_direction',struct('do_stretch',true,'num_samp',1000));
        pos = [td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2)];
        % rotate position by
        R = [cos(-td_temp.target_direction) -sin(-td_temp.target_direction); sin(-td_temp.target_direction) cos(-td_temp.target_direction)];
        newPos = zeros(size(pos));
        for j = 1:length(pos)
            newPos(j,:) = R*(pos(j,:)');
        end
        newPos = zeros(size(pos));
        for j = 1:length(pos)
            newPos(j,:) = R*(pos(j,:)');
        end
        if td(1).perturbation_info(2) < 0
            newPos = -newPos;
        end
        newPos = newPos - repmat(newPos(1,:),size(newPos,1),1);
        plot(newPos(:,1),newPos(:,2),'k-','LineWidth',2);
        
        
        %%%%%%%%%%%%%%%%
        td_temp = trialAverage(td(getTDidx(td,'epoch','AD','target_direction',ad_dir,'range',{'last',num_reaches})), ...
            'target_direction',struct('do_stretch',true,'num_samp',1000));
        pos = [td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2)];
        % rotate position by
        R = [cos(-td_temp.target_direction) -sin(-td_temp.target_direction); sin(-td_temp.target_direction) cos(-td_temp.target_direction)];
        newPos = zeros(size(pos));
        for j = 1:length(pos)
            newPos(j,:) = R*(pos(j,:)');
        end
        newPos = zeros(size(pos));
        for j = 1:length(pos)
            newPos(j,:) = R*(pos(j,:)');
        end
        if td(1).perturbation_info(2) < 0
            newPos = -newPos;
        end
        newPos = newPos - repmat(newPos(1,:),size(newPos,1),1);
        plot(newPos(:,1),newPos(:,2),'b--','LineWidth',2);
        
        
        %%%%%%%%%%%%%%%%%%%
        td_temp = trialAverage(td(getTDidx(td,'epoch','WO','target_direction',ad_dir,'range',{'first',num_reaches})), ...
            'target_direction',struct('do_stretch',true,'num_samp',1000));
        pos = [td_temp.pos(:,1)-pos_offsets(1),td_temp.pos(:,2)-pos_offsets(2)];
        % rotate position by
        R = [cos(-td_temp.target_direction) -sin(-td_temp.target_direction); sin(-td_temp.target_direction) cos(-td_temp.target_direction)];
        newPos = zeros(size(pos));
        for j = 1:length(pos)
            newPos(j,:) = R*(pos(j,:)');
        end
        newPos = zeros(size(pos));
        for j = 1:length(pos)
            newPos(j,:) = R*(pos(j,:)');
        end
        if td(1).perturbation_info(2) < 0
            newPos = -newPos;
        end
        newPos = newPos - repmat(newPos(1,:),size(newPos,1),1);
        plot(newPos(:,1),newPos(:,2),'r-','LineWidth',2);
    end
    legend({'First AD','BL','Last AD','First WO'});
    
    axis('square');
    set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-9 9],'YLim',[-9,9]);
    
    
end
