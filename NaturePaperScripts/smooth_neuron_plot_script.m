%% load data
clear;
close all;
clc;

load('/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-07.mat');

m1_units = [1,3,7,9,28];
units = [11,13,20,24,33];

% filter data
td = removeBadTrials(trial_data,struct('trial_time',[180,280]));
[~,td] = getTDidx(td,'epoch',{'BL'},'result','R');
td = smoothSignals(td,struct('signals',{{'M1_spikes','PMd_spikes'}},'calc_fr',true,'kernel_SD',0.1));
% process data
td1 = trimTD(td,{'idx_target_on',0},{'idx_target_on',50});
td1 = trialAverage(td1,{'target_direction','epoch'},struct('do_stretch',false));
td2 = trimTD(td,{'idx_go_cue',-20},{'idx_go_cue',30});
td2 = trialAverage(td2,{'target_direction','epoch'},struct('do_stretch',false));
td3 = trimTD(td,{'idx_movement_on',0},{'idx_movement_on',50});
td3 = trialAverage(td3,{'target_direction','epoch'},struct('do_stretch',false));

td1 = softNormalize(td1);
td2 = softNormalize(td2);
td3 = softNormalize(td3);

% td1 = subtractConditionMean(td1);
% td2 = subtractConditionMean(td2);
% td3 = subtractConditionMean(td3);
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

for j = 1:length(units)
    array = 'PMd';
    unit = units(j);
    
    subplot1((j-1)*7+1); hold all;
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
    
    subplot1((j-1)*7+2); hold all;
    for i = 1:size(fr2,2)
        plot(fr2(:,i),'Color',c(i,:),'LineWidth',2);
    end
    axis('tight');
    set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
    if j == 1
        title('PMd','FontSize',16);
    end
    
    subplot1((j-1)*7+3); hold all;
    for i = 1:size(fr3,2)
        plot(fr3(:,i),'Color',c(i,:),'LineWidth',2);
    end
    axis('tight');
    set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    array = 'M1';
    unit = m1_units(j);
    
    subplot1((j-1)*7+5); hold all;
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
    
    subplot1((j-1)*7+6); hold all;
    for i = 1:size(fr2,2)
        plot(fr2(:,i),'Color',c(i,:),'LineWidth',2);
    end
    axis('tight');
    set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
    if j == 1
        title('M1','FontSize',16);
    end
    
    subplot1((j-1)*7+7); hold all;
    for i = 1:size(fr3,2)
        plot(fr3(:,i),'Color',c(i,:),'LineWidth',2);
    end
    axis('tight');
    set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
end


%% plot pre/post learning for one target
clear;
clc;

load('/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-07.mat');
% filter data
td = removeBadTrials(trial_data,struct('trial_time',[180,280]));
u = unique([td.target_direction]);
[~,td1] = getTDidx(td,'epoch',{'BL'},'result','R');
[~,td2] = getTDidx(td,'epoch',{'AD'},'range',[0.5 1],'result','R');
td = [td1,td2];

td = smoothSignals(td,struct('signals',{{'M1_spikes','PMd_spikes'}},'calc_fr',true,'kernel_SD',0.1));
% process data
td1 = trimTD(td,{'idx_target_on',0},{'idx_target_on',50});
td1 = trialAverage(td1,{'target_direction','epoch'},struct('do_stretch',false));
td2 = trimTD(td,{'idx_go_cue',-20},{'idx_go_cue',30});
td2 = trialAverage(td2,{'target_direction','epoch'},struct('do_stretch',false));
td3 = trimTD(td,{'idx_movement_on',0},{'idx_movement_on',50});
td3 = trialAverage(td3,{'target_direction','epoch'},struct('do_stretch',false));

% td1 = softNormalize(td1);
% td2 = softNormalize(td2);
% td3 = softNormalize(td3);

% td1 = subtractConditionMean(td1,struct('cond_idx',getTDidx(td1,'epoch','BL')));
% td2 = subtractConditionMean(td2,struct('cond_idx',getTDidx(td2,'epoch','BL')));
% td3 = subtractConditionMean(td3,struct('cond_idx',getTDidx(td3,'epoch','BL')));

% plot stuff
c = phasemap(length(u));

line_style = {'--','--','-','-'};
m1_units = [1,3,7,9,28];
units = [11,13,20,24,33];

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
    
    for j = 1:length(units)
        array = 'PMd';
        unit = units(j);
        
        subplot1((j-1)*7+1); hold all;
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
        
        subplot1((j-1)*7+2); hold all;
        for i = 1:size(fr2,2)
            plot(fr2(:,i),line_style{i},'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
        if j == 1
            title('PMd','FontSize',16);
        end
        
        subplot1((j-1)*7+3); hold all;
        for i = 1:size(fr3,2)
            plot(fr3(:,i),line_style{i},'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        array = 'M1';
        unit = m1_units(j);
        
        subplot1((j-1)*7+5); hold all;
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
        
        subplot1((j-1)*7+6); hold all;
        for i = 1:size(fr2,2)
            plot(fr2(:,i),line_style{i},'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
        if j == 1
            title('M1','FontSize',16);
        end
        
        subplot1((j-1)*7+7); hold all;
        for i = 1:size(fr3,2)
            plot(fr3(:,i),line_style{i},'Color',c(i,:),'LineWidth',2);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[y_min,y_max],'YTick',[]);
    end

