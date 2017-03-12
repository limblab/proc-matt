%% Prep datasets and get paths etc
clear; clc; close all;
dataSummary;
trial_params;
% Session parameters
monkey = 'Chewie';
task = 'CO';
pert = 'FF';
date = '2016-10-07';

units = [20,24,33; ...
    3,9,28];

file = getFileDBidx(filedb, ...
    {'Task',task,'Perturbation',pert,'Monkey',monkey,'Date',date}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});

func_calls = [trial_func_calls, { ...
    {@getTDidx,'epoch',{'BL','AD'}}, ...
    {@smoothSignals,struct('signals',{{'M1_spikes','PMd_spikes'}},'calc_fr',true,'kernel_SD',pn_kernel_SD)}, ...
    {@trimTD,idx_start,idx_end}}];

% load it
fname = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
[trial_data,params] = loadTDfiles(fname, func_calls{:});

%%
[~,td] = getTDidx(trial_data,'epoch','BL');
td1 = trimTD(td,{'idx_target_on',0},{'idx_target_on',60});
td1 = trialAverage(td1,{'target_direction','epoch'},struct('do_stretch',false));
td2 = trimTD(td,{'idx_go_cue',-25},{'idx_go_cue',25});
td2 = trialAverage(td2,{'target_direction','epoch'},struct('do_stretch',false));
td3 = trimTD(td,{'idx_movement_on',0},{'idx_movement_on',50});
td3 = trialAverage(td3,{'target_direction','epoch'},struct('do_stretch',false));

% plot stuff
close all;
u = unique([td.target_direction]);
c = phasemap(length(u));

figure;
subplot1(length(units),7);

for j = 1:length(units)
    subplot1((j-1)*7+4);
    set(gca,'Visible','off');
end

arrays = {'PMd','M1'};
for a = 1:length(arrays)
    array = arrays{a};
    for j = 1:length(units)
        unit = units(a,j);
        
        subplot1((j-1)*7+1+4*(a-1)); hold all;
        fr1 = cell2mat(cellfun(@(x) x(:,unit),{td1.([array '_spikes'])},'Uni',0));
        fr2 = cell2mat(cellfun(@(x) x(:,unit),{td2.([array '_spikes'])},'Uni',0));
        fr3 = cell2mat(cellfun(@(x) x(:,unit),{td3.([array '_spikes'])},'Uni',0));
        
        y_min = min([min(min(fr1)), min(min(fr2)), min(min(fr3))]);
        y_max = max([max(max(fr1)), max(max(fr2)), max(max(fr3))]);
        
        for i = 1:size(fr1,2)
            plot(fr1(:,i),'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
        
        subplot1((j-1)*7+2+4*(a-1)); hold all;
        for i = 1:size(fr2,2)
            plot(fr2(:,i),'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
        if j == 1
            title('PMd','FontSize',16);
        end
        
        subplot1((j-1)*7+3+4*(a-1)); hold all;
        for i = 1:size(fr3,2)
            plot(fr3(:,i),'Color',c(i,:),'LineWidth',2);
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
td1 = trimTD(td,{'idx_target_on',0},{'idx_target_on',60});
td1 = trialAverage(td1,{'target_direction','epoch'},struct('do_stretch',false));
td2 = trimTD(td,{'idx_go_cue',-25},{'idx_go_cue',25});
td2 = trialAverage(td2,{'target_direction','epoch'},struct('do_stretch',false));
td3 = trimTD(td,{'idx_movement_on',0},{'idx_movement_on',50});
td3 = trialAverage(td3,{'target_direction','epoch'},struct('do_stretch',false));

line_style = {'--','--','-','-'};

figure;
subplot1(length(units),7);

which_target = [2,6];
c = c(which_target([1 2 1 2]),:);

[~,td1] = getTDidx(td1,'target_direction',u(which_target));
[~,td2] = getTDidx(td2,'target_direction',u(which_target));
[~,td3] = getTDidx(td3,'target_direction',u(which_target));

for j = 1:length(units)
    subplot1((j-1)*7+4);
    set(gca,'Visible','off');
end
arrays = {'PMd','M1'};
for a = 1:length(arrays)
    array = arrays{a};
    for j = 1:size(units,2)
        
        unit = units(a,j);
        
        subplot1((j-1)*7+1+4*(a-1)); hold all;
        fr1 = cell2mat(cellfun(@(x) x(:,unit),{td1.([array '_spikes'])},'Uni',0));
        fr2 = cell2mat(cellfun(@(x) x(:,unit),{td2.([array '_spikes'])},'Uni',0));
        fr3 = cell2mat(cellfun(@(x) x(:,unit),{td3.([array '_spikes'])},'Uni',0));
        
        y_min = min([min(min(fr1)), min(min(fr2)), min(min(fr3))]);
        y_max = max([max(max(fr1)), max(max(fr2)), max(max(fr3))]);
        
        for i = 1:size(fr1,2)
            plot(fr1(:,i),line_style{i},'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
        
        subplot1((j-1)*7+2+4*(a-1)); hold all;
        for i = 1:size(fr2,2)
            plot(fr2(:,i),line_style{i},'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
        if j == 1
            title('PMd','FontSize',16);
        end
        
        subplot1((j-1)*7+3+4*(a-1)); hold all;
        for i = 1:size(fr3,2)
            plot(fr3(:,i),line_style{i},'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
    end
end

