%% Prep datasets and get paths etc
clear; clc; close all;
dataSummary;

% Session parameters
monkey = 'Chewie';
task = 'CO';
pert = 'FF';
date = '2016-10-07';

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
td1 = trimTD(td,{'idx_target_on',-20},{'idx_target_on',60});
td1 = trialAverage(td1,{'target_direction','epoch'},struct('do_stretch',false));
td2 = trimTD(td,{'idx_go_cue',-20},{'idx_go_cue',80});
td2 = trialAverage(td2,{'target_direction','epoch'},struct('do_stretch',false));

% plot stuff
close all;
u = unique([td.target_direction]);
c = phasemap(length(u));

figure;
subplot1(length(units),5);

for j = 1:length(units)
    subplot1((j-1)*5+3);
    set(gca,'Visible','off');
end

arrays = {'PMd','M1'};
for a = 1:length(arrays)
    array = arrays{a};
    for j = 1:length(units)
        unit = units(a,j);
        
        subplot1((j-1)*5+1+3*(a-1)); hold all;
        fr1 = cell2mat(cellfun(@(x) x(:,unit),{td1.([array '_spikes'])},'Uni',0));
        fr2 = cell2mat(cellfun(@(x) x(:,unit),{td2.([array '_spikes'])},'Uni',0));
        
        y_min = min([min(min(fr1)), min(min(fr2))]);
        y_max = max([max(max(fr1)), max(max(fr2))]);
        
        for i = 1:size(fr1,2)
            plot(fr1(:,i),'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
        
                if j == 1
            title(array,'FontSize',16);
                end
        
        subplot1((j-1)*5+2+3*(a-1)); hold all;
        for i = 1:size(fr2,2)
            plot(fr2(:,i),'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
        
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
td1 = trimTD(td,{'idx_target_on',-20},{'idx_target_on',60});
td1 = trialAverage(td1,{'target_direction','epoch'},struct('do_stretch',false));
td2 = trimTD(td,{'idx_go_cue',-20},{'idx_go_cue',80});
td2 = trialAverage(td2,{'target_direction','epoch'},struct('do_stretch',false));

line_style = {'--','--','-','-'};

figure;
subplot1(length(units),5);

which_target = [2,6];
c = c(which_target([1 2 1 2]),:);

[~,td1] = getTDidx(td1,'target_direction',u(which_target));
[~,td2] = getTDidx(td2,'target_direction',u(which_target));

for j = 1:length(units)
    subplot1((j-1)*5+3);
    set(gca,'Visible','off');
end
arrays = {'PMd','M1'};
for a = 1:length(arrays)
    array = arrays{a};
    for j = 1:size(units,2)
        
        unit = units(a,j);
        
        subplot1((j-1)*5+1+3*(a-1)); hold all;
        fr1 = cell2mat(cellfun(@(x) x(:,unit),{td1.([array '_spikes'])},'Uni',0));
        fr2 = cell2mat(cellfun(@(x) x(:,unit),{td2.([array '_spikes'])},'Uni',0));
        
        y_min = min([min(min(fr1)), min(min(fr2))]);
        y_max = max([max(max(fr1)), max(max(fr2))]);
        
        for i = 1:size(fr1,2)
            plot(fr1(:,i),line_style{i},'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
        
        subplot1((j-1)*5+2+3*(a-1)); hold all;
        for i = 1:size(fr2,2)
            plot(fr2(:,i),line_style{i},'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
        if j == 1
            title(array,'FontSize',16);
        end
    end
end

