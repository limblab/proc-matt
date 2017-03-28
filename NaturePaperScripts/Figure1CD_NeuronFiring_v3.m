% THIS VERSION ONLY HAS ONE WINDOW CENTERED ON GO CUE
%% Prep datasets and get paths etc
clear; clc; close all;
dataSummary;

% Session parameters
monkey = 'Chewie';
task = 'CO';
pert = 'VR';
date = '2016-10-06';

units = [20,24,33; ...
    3,9,28];

badneuron_params = struct( ...
    'min_fr',0.3, ...
    'do_shunt_check',1, ...
    'use_trials',{{'epoch','BL'}});

file = getFileDBidx(filedb, ...
    {'Task',task,'Perturbation',pert,'Monkey',monkey,'Date',date}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});

idx_start = {'idx_target_on', -20};
idx_end   = {'idx_trial_end', -20};

func_calls = { ...
    {@getTDidx,'result','r'}, ...
    {@removeBadNeurons,badneuron_params}, ...
    {@getTDidx,'epoch',{'BL','AD'}}, ...
    {@smoothSignals,struct('signals',{{'M1_spikes','PMd_spikes'}},'calc_fr',true,'kernel_SD',0.1)}, ...
    {@trimTD,idx_start,idx_end}};

% load it
fname = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
[trial_data,params] = loadTDfiles(fname, func_calls{:});

%%
[~,td] = getTDidx(trial_data,'epoch','BL');
td = trimTD(td,{'idx_go_cue',-70},{'idx_go_cue',80});
td = trialAverage(td,{'target_direction','epoch'},struct('do_stretch',false));

% plot stuff
close all;
u = unique([td.target_direction]);
c = phasemap(length(u));

figure;
subplot1(length(units),2);

arrays = {'PMd','M1'};
for a = 1:length(arrays)
    array = arrays{a};
    for j = 1:length(units)
        unit = units(a,j);
        
        subplot1(2*(j-1)+a); hold all;
        fr = cell2mat(cellfun(@(x) x(:,unit),{td.([array '_spikes'])},'Uni',0));
        
        y_min = min(min(fr));
        y_max = max(max(fr));
        
        for i = 1:size(fr,2)
            plot(fr(:,i),'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
        
        if j == 1
            title(array,'FontSize',16);
        end
        
        
    end
end

%%
td = trial_data;

% filter data
u = unique([td.target_direction]);
[~,td1] = getTDidx(td,'epoch',{'BL'},'result','R');
[~,td2] = getTDidx(td,'epoch',{'AD'},'range',[0.5 1],'result','R');
td = [td1,td2];

% process data
td = trimTD(td,{'idx_go_cue',-70},{'idx_go_cue',80});
td = trialAverage(td,{'target_direction','epoch'},struct('do_stretch',false));

line_style = {'--','--','-','-'};

figure;
subplot1(length(units),2);

which_target = [2,6];
c = c(which_target([1 2 1 2]),:);

[~,td] = getTDidx(td,'target_direction',u(which_target));

arrays = {'PMd','M1'};
for a = 1:length(arrays)
    array = arrays{a};
    for j = 1:size(units,2)
        
        unit = units(a,j);
        
        subplot1(2*(j-1)+a); hold all;
        fr = cell2mat(cellfun(@(x) x(:,unit),{td.([array '_spikes'])},'Uni',0));
        
        y_min = min(min(fr));
        y_max = max(max(fr));
        
        for i = 1:size(fr,2)
            plot(fr(:,i),line_style{i},'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
        
        if j == 1
            title(array,'FontSize',16);
        end
    end
end

