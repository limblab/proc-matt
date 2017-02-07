%% CELL CLASSES AND SCATTER PLOT
clear;
close all;
clc;

dataSummary;

use_model = 'muscle';
reassign_others = false;

alpha = 0.05;
comp_blocks = [1 2 3];

min_fr = 5;
min_r2 = 0.5;
min_weights = -Inf;
max_weights = Inf;
num_plot_neurons = 'all';
tune_real_neurons = false;

classColors = {'k','b','r','m','g','c'};
% root_dir = 'F:\m1_cf_JNS_results\biomech_sim_results\';
dataDir = fullfile(rootDir,resultsDir,'biomech_model');

filenames = { ...
%     'Chewie_CO_FF_2013-10-22', ...
%     'Chewie_CO_FF_2013-10-23', ...
%     'Chewie_CO_FF_2013-10-31', ...
%     'Chewie_CO_FF_2013-11-01', ...
%     'Chewie_CO_FF_2013-12-03', ...
%         'Chewie_CO_FF_2013-12-04', ...
%     'Chewie_CO_FF_2015-07-01', ...
    'Chewie_CO_FF_2015-07-03', ...
%     'Mihili_CO_FF_2014-02-03', ...
%     'Mihili_CO_FF_2014-02-17', ...
%     'Mihili_CO_FF_2014-02-18', ...
%     'Mihili_CO_FF_2014-03-07', ...
%     'Mihili_CO_FF_2015-06-11', ...
%     'Mihili_CO_FF_2015-06-17', ...
%     'Chewie_CO_FF_2015-06-29', ...
%     'Chewie_CO_FF_2015-06-30', ...
%     'Chewie_CO_FF_2015-07-06', ...
%         'Chewie_CO_FF_2015-07-07', ...
%             'Chewie_CO_FF_2015-07-08', ...
%     'Mihili_CO_FF_2015-06-10', ...
%     'Mihili_CO_FF_2015-06-15', ...
%     'Mihili_CO_FF_2015-06-16', ...
    };

[dpd_ad, dpd_wo, dpd_adwo, ddom_ad, cell_classes, bl_avg, bl_dom, bl_r2, bl_cb,bl_pd,tuned_cells,cell_tuning_idx,cell_tcs,all_weights] = deal([]);
num_cells = zeros(1,length(filenames));
for iFile = 1:length(filenames)
    load(fullfile(dataDir,[filenames{iFile} '_tuning.mat']),'tc_data','params');
    load(fullfile(dataDir,[filenames{iFile} '_neurons.mat']),'neural_tcs');
    
    num_neurons = size(tc_data(1).muscle.tc,1);
    
    % get index describing tuning of each cell
    % NOTE: DOESN'T WORK WITH KIN YET
    tc = neural_tcs.(use_model);
    tc_index = zeros(size(tc,1),2);
    for i = 1:size(tc,1)
        if tc(i,2) + tc(i,4) ~= 0 && tc(i,3) + tc(i,5) == 0 % it's a flexor cell
            % 0 means equal weight to both joints
            % positive means more shoulder
            % negative means more elbow
            tc_index(i,1) = 1;
            tc_index(i,2) = (abs(tc(i,2)) - abs(tc(i,4)))/(abs(tc(i,2)) + abs(tc(i,4)));
        elseif tc(i,3) + tc(i,5) ~= 0 && tc(i,2) + tc(i,4) == 0 % it's an extensor cell
            tc_index(i,1) = 2;
            tc_index(i,2) = (abs(tc(i,3)) - abs(tc(i,5)))/(abs(tc(i,3)) + abs(tc(i,5)));
            
        elseif tc(i,2) + tc(i,3) ~= 0 && tc(i,4) + tc(i,5) == 0 % it's a shoulder
            tc_index(i,1) = 3;
            tc_index(i,2) = (abs(tc(i,2)) - abs(tc(i,3)))/(abs(tc(i,2)) + abs(tc(i,3)));
        elseif tc(i,4) + tc(i,5) ~= 0 && tc(i,2) + tc(i,3) == 0 % it's an elbow cell
            tc_index(i,1) = 4;
            tc_index(i,2) = (abs(tc(i,4)) - abs(tc(i,5)))/(abs(tc(i,4)) + abs(tc(i,5)));
        else % we didn't do synergies or joints
            % OPTION 1
            % 0 means equal weight to flexion and extension
            % negative means more flexion
            % positive means more extension
            % tc_index(i,2) = ((abs(tc(i,2))+abs(tc(i,4))) - (abs(tc(i,3))+abs(tc(i,5))))/((abs(tc(i,2))+abs(tc(i,4))) + (abs(tc(i,3))+abs(tc(i,5))));
            
            % OPTION 2
            % 0 means equal elbow and shoulder
            % negative means more shoulder
            % positive means more elbow
            tc_index(i,2) = ((abs(tc(i,2))+abs(tc(i,3))) - (abs(tc(i,4))+abs(tc(i,5))))/((abs(tc(i,2))+abs(tc(i,3))) + (abs(tc(i,4))+abs(tc(i,5))));
        end
    end
    
    if ~tune_real_neurons
        num_cells(iFile) = sum(mean(tc_data(1).muscle.boot_bos,2) > min_fr & sum(abs(neural_tcs.muscle(:,2:end)),2) > min_weights & sum(abs(neural_tcs.muscle(:,2:end)),2) < max_weights & any(neural_tcs.muscle(:,2:end) > 0,2));
        all_weights = [all_weights; sum(abs(neural_tcs.muscle(:,2:end)),2)];
    end
    
    %     mean(tc_data(1).(use_model).rs,2) > 0.5 & ...
    is_tuned = all(prctile(tc_data(comp_blocks(1)).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
        all(prctile(tc_data(comp_blocks(2)).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
        all(prctile(tc_data(comp_blocks(3)).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
        angleDiff(tc_data(comp_blocks(1)).(use_model).cb{3}(:,1),tc_data(comp_blocks(1)).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
        angleDiff(tc_data(comp_blocks(2)).(use_model).cb{3}(:,1),tc_data(comp_blocks(2)).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
        angleDiff(tc_data(comp_blocks(3)).(use_model).cb{3}(:,1),tc_data(comp_blocks(3)).(use_model).cb{3}(:,2),true,false) < 40*pi/180;
    
    % get some facts about the tuned cells
    bl_avg = [bl_avg; tc_data(1).(use_model).tc(:,1)]; % mean
    bl_dom = [bl_dom; tc_data(1).(use_model).tc(:,2)]; % modulation depth
    bl_r2 = [bl_r2; mean(tc_data(1).(use_model).rs,2)]; % r-squared
    bl_cb = [bl_cb; angleDiff(tc_data(1).(use_model).cb{3}(:,1),tc_data(1).(use_model).cb{3}(:,2),true,false)];
    bl_pd = [bl_pd; tc_data(1).(use_model).tc(:,3)];
    
    % get PDs
    bl = tc_data(comp_blocks(1)).(use_model).tc(:,3);
    ad = tc_data(comp_blocks(2)).(use_model).tc(:,3);
    wo = tc_data(comp_blocks(3)).(use_model).tc(:,3);
    
    temp = angleDiff(bl,ad,true,true);
    
        if mean(temp) < 0
            disp('FLIPPING');
            temp = -temp;
        end
%     temp = sign(trial_data(1).perturbation_info(2))*temp;
    
    dpd_ad = [dpd_ad; temp];
    dpd_wo = [dpd_wo; angleDiff(bl,wo,true,true)];
    dpd_adwo = [dpd_adwo; angleDiff(ad,wo,true,true)];
    
    % get DOM
    bl = tc_data(comp_blocks(1)).(use_model).tc(:,2);
    ad = tc_data(comp_blocks(2)).(use_model).tc(:,2);
    wo = tc_data(comp_blocks(3)).(use_model).tc(:,2);
    ddom_ad = [ddom_ad; ad-bl];
    
    % classify
    all_perms = nchoosek(comp_blocks,2);
    
    is_diff = zeros(num_neurons,size(all_perms,1));
    for j = 1:size(all_perms,1)
        for i = 1:size(tc_data(1).(use_model).boot_pds,1)
            cb = prctile(angleDiff(tc_data(all_perms(j,1)).(use_model).boot_pds(i,:),tc_data(all_perms(j,2)).(use_model).boot_pds(i,:),true,true),100*[alpha/2,1-alpha/2],2);
            %             cb = prctile(tc_data(all_perms(j,2)).(use_model).boot_pds(i,:)-tc_data(all_perms(j,1)).(use_model).boot_pds(i,:),100*[alpha/2,1-alpha/2],2);
            
            if isempty(range_intersection([0 0],cb)) || (abs(temp(i)) > 150*pi/180 && j == 1) || (abs(temp(i)) > 150*pi/180 && j == 3)
                is_diff(i,j) = 1;
            end
        end
    end
    
    cc = zeros(size(is_diff,1),1);
    for i = 1:size(is_diff,1)
        % 2 dynamic: 1 0 1
        % 1 kinematic: 0 0 0
        % 3 memory I: 1 1 0
        % 4 memory II: 0 1 1
        % 5 other: 1 1 1
        if all(is_diff(i,:) == [0 0 0])
            cc(i) = 1;
        elseif all(is_diff(i,:) == [1 0 1])
            cc(i) = 2;
        elseif all(is_diff(i,:) == [1 1 0]) || all(is_diff(i,:) == [1 0 0])
            cc(i) = 3;
        elseif all(is_diff(i,:) == [0 1 1]) || all(is_diff(i,:) == [0 0 1]) || all(is_diff(i,:) == [0 1 0])
            cc(i) = 4;
        elseif all(is_diff(i,:) == [1 1 1])
            cc(i) = 5;
        else
            disp('fuck');
            cc(i) = 6;
        end
    end
    
    cell_classes = [cell_classes; cc];
    tuned_cells = [tuned_cells; is_tuned];
    cell_tuning_idx = [cell_tuning_idx; tc_index];
    cell_tcs = [cell_tcs; tc];
end

if reassign_others
    mem_ind = abs(dpd_wo) ./ min( abs(dpd_ad) , abs(dpd_adwo) );
    idx = cell_classes == 5 & mem_ind >= 1;
    cell_classes(idx) = 3;
    idx = cell_classes == 5 & mem_ind < 1;
    cell_classes(idx) = 2;
end

% tuned_cells = ~tuned_cells;

tuned_cells = logical(tuned_cells);

% if ~tune_real_neurons
%     tuned_cells = tuned_cells & ...
%         ~(cell_tcs(:,2) < 0 & cell_tcs(:,4) > 0 & cell_tcs(:,6) < 0) & ...
%         ~(cell_tcs(:,3) < 0 & cell_tcs(:,5) > 0 & cell_tcs(:,7) > 0);
% end

bl_pd = bl_pd(tuned_cells);
bl_avg = bl_avg(tuned_cells);
bl_dom = bl_dom(tuned_cells);
bl_r2 = bl_r2(tuned_cells);
bl_cb = bl_cb(tuned_cells);
dpd_ad = dpd_ad(tuned_cells);
dpd_wo = dpd_wo(tuned_cells);
dpd_adwo = dpd_adwo(tuned_cells);
ddom_ad = ddom_ad(tuned_cells);
cell_classes = cell_classes(tuned_cells);

cell_tuning_joint = cell_tuning_idx(tuned_cells,1);
cell_tuning_idx = cell_tuning_idx(tuned_cells,2);
cell_tcs = cell_tcs(tuned_cells,:);

binsize = 5;
figure('Position',[300 50 950 950]);
subplot1(2,2,'Gap',[0 0]);
subplot1(3);
hold all;

if strcmpi(num_plot_neurons,'all')
    temp_dpd = 1:length(dpd_ad);
else
    temp_dpd = randi(length(dpd_ad),1,num_plot_neurons);
end

for i  = temp_dpd
    plot(dpd_ad(i).*(180/pi),dpd_wo(i).*(180/pi),'d','LineWidth',2,'Color',classColors{cell_classes(i)});
end
plot([-180,180],[0 0],'k--','LineWidth',1);
plot([0 0],[-180,180],'k--','LineWidth',1);
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-180,180],'YLim',[-180,180]);

xlabel('dPD Base to Force','FontSize',14);
ylabel('dPD Base to Wash','FontSize',14);

subplot1(1);
hold all;
hist(dpd_ad.*180/pi,-180:binsize:180);
axis('tight');
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-180,180],'XTick',[],'YTick',[]);

subplot1(4);
hold all;
hist(dpd_wo.*180/pi,-180:binsize:180);
axis('tight');
set(gca,'CameraUpVector',[1,0,0],'Xdir','reverse');
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-180,180],'XTick',[],'YTick',[]);

subplot1(2);
set(gca,'Box','off','Visible','off');

figure;
[k,x] = hist(cell_classes,1:5);
bar(x,100*k/sum(k),1)
set(gca,'Box','off','TickDir','out','FontSize',14,'XTick',1:5,'XTickLabel',{'Kin','Dyn','MemI','MemII','Other'},'XLim',[0 6]);
ylabel('Percent of Cells','FontSize',14);
title([use_model '-based neurons'],'FontSize',14);

%%
figure;
subplot(2,3,1); plot(dpd_ad, cell_tcs(:,2), '.'); title('shoulder flex');
subplot(2,3,2); plot(dpd_ad, cell_tcs(:,4), '.'); title('elbow flex');
subplot(2,3,3); plot(dpd_ad, cell_tcs(:,6), '.'); title('biart flex');
subplot(2,3,4); plot(dpd_ad, cell_tcs(:,3), '.'); title('shoulder ext');
subplot(2,3,5); plot(dpd_ad, cell_tcs(:,5), '.'); title('elbow ext');
subplot(2,3,6); plot(dpd_ad, cell_tcs(:,7), '.'); title('biart ext');

%%
% close all;
p = anovan(dpd_ad,[cell_tcs(:,2:end), ...
    abs(cell_tcs(:,3)-cell_tcs(:,2)), ...
    abs(cell_tcs(:,5)-cell_tcs(:,4)), ...
    abs(cell_tcs(:,4)-cell_tcs(:,2)), ...
    abs(cell_tcs(:,5)-cell_tcs(:,3)), ...
    ddom_ad, ...
    bl_avg, ...
    bl_dom, ...
    bl_r2],'continuous',1:12);

figure;
subplot(3,4,1); plot(dpd_ad, cell_tcs(:,2), '.'); title(['shoulder flex; p=' num2str(p(1))]);
subplot(3,4,2); plot(dpd_ad, cell_tcs(:,3), '.'); title(['shoulder ext; p=' num2str(p(2))]);
subplot(3,4,3); plot(dpd_ad, cell_tcs(:,4), '.'); title(['elbow flex; p=' num2str(p(3))]);
subplot(3,4,4); plot(dpd_ad, cell_tcs(:,5), '.'); title(['elbow ext; p=' num2str(p(4))]);

subplot(3,4,5); plot(dpd_ad, abs(cell_tcs(:,3)-cell_tcs(:,2)), '.'); title(['shoulder flex/ext; p=' num2str(p(5))]);
subplot(3,4,6); plot(dpd_ad, abs(cell_tcs(:,5)-cell_tcs(:,4)), '.'); title(['elbow flex/ext; p=' num2str(p(6))]);
subplot(3,4,7); plot(dpd_ad, abs(cell_tcs(:,4)-cell_tcs(:,2)), '.'); title(['elbow/shoulder flex; p=' num2str(p(7))]);
subplot(3,4,8); plot(dpd_ad, abs(cell_tcs(:,5)-cell_tcs(:,3)), '.'); title(['elbow/shoulder ext; p=' num2str(p(8))]);

subplot(3,4,9); plot(dpd_ad, ddom_ad, '.'); title(['change in dom; p=' num2str(p(9))]);
subplot(3,4,10); plot(dpd_ad, bl_avg, '.'); title(['baseline fr; p=' num2str(p(10))]);
subplot(3,4,11); plot(dpd_ad, bl_dom, '.'); title(['baseline dom; p=' num2str(p(11))]);
subplot(3,4,12); plot(dpd_ad, bl_r2, '.'); title(['baseline r2; p=' num2str(p(12))]);


%% OVER TIME

use_model = 'muscle';
comp_blocks = [1 4 7];
errScale = 2;

% root_dir = 'F:\trial_data_files\biomech_sim_results\';

filenames = { ...
    'Chewie_CO_FF_2013-10-22', ...
    'Chewie_CO_FF_2013-10-23', ...
    'Chewie_CO_FF_2013-10-31', ...
    'Chewie_CO_FF_2013-11-01', ...
    'Chewie_CO_FF_2013-12-03', ...
        'Chewie_CO_FF_2013-12-04', ...
    'Chewie_CO_FF_2015-07-01', ...
    'Chewie_CO_FF_2015-07-03', ...
    'Chewie_CO_FF_2015-06-29', ...
    'Chewie_CO_FF_2015-06-30', ...
    'Chewie_CO_FF_2015-07-06', ...
            'Chewie_CO_FF_2015-07-07', ...
                    'Chewie_CO_FF_2015-07-08', ...
    };


[dpd,dpd_err,all_bl] = deal([]);
for iFile = 1:length(filenames)
    load([dataDir filenames{iFile} '_tuning.mat'],'tc_data','params');
    
    %     mean(tc_data(1).(use_model).rs,2) > 0.5 & ...
    is_tuned = all(prctile(tc_data(comp_blocks(1)).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
        all(prctile(tc_data(comp_blocks(2)).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
        all(prctile(tc_data(comp_blocks(3)).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
        angleDiff(tc_data(comp_blocks(1)).(use_model).cb{3}(:,1),tc_data(comp_blocks(1)).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
        angleDiff(tc_data(comp_blocks(2)).(use_model).cb{3}(:,1),tc_data(comp_blocks(2)).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
        angleDiff(tc_data(comp_blocks(3)).(use_model).cb{3}(:,1),tc_data(comp_blocks(3)).(use_model).cb{3}(:,2),true,false) < 40*pi/180;
    
    % get PDs
    bl = tc_data(comp_blocks(1)).(use_model).tc(:,3);
    
    file_dpd = zeros(sum(is_tuned),7);
    file_dpd_err = zeros(sum(is_tuned),7);
    for i = 1:7
        temp = angleDiff(bl,tc_data(i).(use_model).tc(:,3),true,true);
        
        if mean(temp) < 0
            temp = -temp;
        end
        file_dpd(:,i) = temp(is_tuned);
        
        temp = angleDiff(tc_data(i).(use_model).cb{3}(:,1),tc_data(i).(use_model).cb{3}(:,2),true,false)./errScale;
        file_dpd_err(:,i) = temp(is_tuned);
    end
    dpd = [dpd; file_dpd];
    dpd_err = [dpd_err; file_dpd_err];
    all_bl = [all_bl; bl];
end

m = circular_mean(dpd,[],1)*180/pi;
s = mean(dpd_err*180/pi);
% s = circular_std(dpd)*180/pi;

figure('Position',[200 200 1280 800]); hold all;
plot(1:7,m,'bo','LineWidth',2);
plot([1:7;1:7],[m-s; m+s],'b-','LineWidth',2);


filenames = { ...
    'Mihili_CO_FF_2014-02-03', ...
    'Mihili_CO_FF_2014-02-17', ...
    'Mihili_CO_FF_2014-02-18', ...
    'Mihili_CO_FF_2014-03-07', ...
    'Mihili_CO_FF_2015-06-11', ...
    'Mihili_CO_FF_2015-06-17', ...
    'Mihili_CO_FF_2015-06-10', ...
    'Mihili_CO_FF_2015-06-15', ...
    'Mihili_CO_FF_2015-06-16', ...
    };


[dpd,dpd_err,all_bl] = deal([]);
for iFile = 1:length(filenames)
    load([dataDir filenames{iFile} '_tuning.mat'],'tc_data','params');
    
    %     mean(tc_data(1).(use_model).rs,2) > 0.5 & ...
    is_tuned = all(prctile(tc_data(comp_blocks(1)).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
        all(prctile(tc_data(comp_blocks(2)).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
        all(prctile(tc_data(comp_blocks(3)).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
        angleDiff(tc_data(comp_blocks(1)).(use_model).cb{3}(:,1),tc_data(comp_blocks(1)).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
        angleDiff(tc_data(comp_blocks(2)).(use_model).cb{3}(:,1),tc_data(comp_blocks(2)).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
        angleDiff(tc_data(comp_blocks(3)).(use_model).cb{3}(:,1),tc_data(comp_blocks(3)).(use_model).cb{3}(:,2),true,false) < 40*pi/180;
    
    % get PDs
    bl = tc_data(comp_blocks(1)).(use_model).tc(:,3);
    
    file_dpd = zeros(sum(is_tuned),7);
    file_dpd_err = zeros(sum(is_tuned),7);
    for i = 1:7
        temp = angleDiff(bl,tc_data(i).(use_model).tc(:,3),true,true);
        
        if mean(temp) < 0
            temp = -temp;
        end
        file_dpd(:,i) = temp(is_tuned);
        
        temp = angleDiff(tc_data(i).(use_model).cb{3}(:,1),tc_data(i).(use_model).cb{3}(:,2),true,false)./errScale;
        file_dpd_err(:,i) = temp(is_tuned);
    end
    dpd = [dpd; file_dpd];
    dpd_err = [dpd_err; file_dpd_err];
    all_bl = [all_bl; bl];
end

plotMin = -40;
plotMax = 110;

m = circular_mean(dpd,[],1)*180/pi;
s = mean(dpd_err*180/pi);
% s = circular_std(dpd)*180/pi;

plot((1:7)+0.1,m,'ro','LineWidth',2);
plot([1:7;1:7]+0.1,[m-s; m+s],'r-','LineWidth',2);


% Now add some extra stuff and configure the plot
plot([0 7+1],[0 0],'LineWidth',1,'Color',[0.6 0.6 0.6]);
plot([1.5 1.5],[plotMin plotMax],'k--','LineWidth',1);
plot([4.5 4.5],[plotMin plotMax],'k--','LineWidth',1);

axis([0.3 7+0.3 plotMin plotMax]);

set(gca,'Box','off','TickDir','out','FontSize',14);



%% CLASSES IN DIFFERENT BLOCKS

use_model = 'muscle';
reassign_others = false;

alpha = 0.05;
class_blocks = {[1 2 5],[1 3 6],[1 4 7]};
% num_neurons = 200;

classColors = {'k','b','r','m','g','c'};

% root_dir = 'F:\trial_data_files\biomech_sim_results\';

filenames = { ...
    'Chewie_CO_FF_2013-10-22', ...
    'Chewie_CO_FF_2013-10-23', ...
    'Chewie_CO_FF_2013-10-31', ...
    'Chewie_CO_FF_2013-11-01', ...
    'Chewie_CO_FF_2013-12-03', ...
        'Chewie_CO_FF_2013-12-04', ...
    'Chewie_CO_FF_2015-07-01', ...
    'Chewie_CO_FF_2015-07-03', ...
    'Mihili_CO_FF_2014-02-03', ...
    'Mihili_CO_FF_2014-02-17', ...
    'Mihili_CO_FF_2014-02-18', ...
    'Mihili_CO_FF_2014-03-07', ...
    'Mihili_CO_FF_2015-06-11', ...
    'Mihili_CO_FF_2015-06-17', ...
    'Chewie_CO_FF_2015-06-29', ...
    'Chewie_CO_FF_2015-06-30', ...
    'Chewie_CO_FF_2015-07-06', ...
            'Chewie_CO_FF_2015-07-07', ...
                    'Chewie_CO_FF_2015-07-08', ...
    'Mihili_CO_FF_2015-06-10', ...
    'Mihili_CO_FF_2015-06-15', ...
    'Mihili_CO_FF_2015-06-16', ...
    };

figure; hold all;
for iBlock = 1:3
    comp_blocks = class_blocks{iBlock};
    [dpd_ad, dpd_wo, dpd_adwo, cell_classes,tuned_cells] = deal([]);
    s_classes = cell(1,length(filenames));
    for iFile = 1:length(filenames)
        load([dataDir filenames{iFile} '_tuning.mat'],'tc_data','neural_tcs','params');
        
        num_neurons = size(tc_data(1).muscle.tc,1);
        
        is_tuned = all(prctile(tc_data(comp_blocks(1)).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
            all(prctile(tc_data(comp_blocks(2)).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
            all(prctile(tc_data(comp_blocks(3)).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
            angleDiff(tc_data(comp_blocks(1)).(use_model).cb{3}(:,1),tc_data(comp_blocks(1)).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
            angleDiff(tc_data(comp_blocks(2)).(use_model).cb{3}(:,1),tc_data(comp_blocks(2)).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
            angleDiff(tc_data(comp_blocks(3)).(use_model).cb{3}(:,1),tc_data(comp_blocks(3)).(use_model).cb{3}(:,2),true,false) < 40*pi/180;
        
        % get PDs
        bl = tc_data(comp_blocks(1)).(use_model).tc(:,3);
        ad = tc_data(comp_blocks(2)).(use_model).tc(:,3);
        wo = tc_data(comp_blocks(3)).(use_model).tc(:,3);
        
        temp = angleDiff(bl,ad,true,true);
        
        if mean(temp) < 0
            temp = -temp;
        end
        
        dpd_ad = [dpd_ad; temp];
        dpd_wo = [dpd_wo; angleDiff(bl,wo,true,true)];
        dpd_adwo = [dpd_adwo; angleDiff(ad,wo,true,true)];
        
        % classify
        all_perms = nchoosek(comp_blocks,2);
        
        is_diff = zeros(num_neurons,size(all_perms,1));
        for j = 1:size(all_perms,1)
            for i = 1:size(tc_data(1).(use_model).boot_pds,1)
                cb = prctile(angleDiff(tc_data(all_perms(j,1)).(use_model).boot_pds(i,:),tc_data(all_perms(j,2)).(use_model).boot_pds(i,:),true,true),100*[alpha/2,1-alpha/2],2);
                %cb = prctile(tc_data(all_perms(j,2)).(use_model).boot_pds(i,:)-tc_data(all_perms(j,1)).(use_model).boot_pds(i,:),100*[alpha/2,1-alpha/2],2);
                
                if isempty(range_intersection([0 0],cb)) || (abs(temp(i)) > 150*pi/180 && j == 1) || (abs(temp(i)) > 150*pi/180 && j == 3)
                    is_diff(i,j) = 1;
                end
            end
        end
        
        cc = zeros(size(is_diff,1),1);
        for i = 1:size(is_diff,1)
            % 2 dynamic: 1 0 1
            % 1 kinematic: 0 0 0
            % 3 memory I: 1 1 0
            % 4 memory II: 0 1 1
            % 5 other: 1 1 1
            if all(is_diff(i,:) == [0 0 0])
                cc(i) = 1;
            elseif all(is_diff(i,:) == [1 0 1])
                cc(i) = 2;
            elseif all(is_diff(i,:) == [1 1 0]) || all(is_diff(i,:) == [1 0 0])
                cc(i) = 3;
            elseif all(is_diff(i,:) == [0 1 1]) || all(is_diff(i,:) == [0 0 1]) || all(is_diff(i,:) == [0 1 0])
                cc(i) = 4;
            elseif all(is_diff(i,:) == [1 1 1])
                cc(i) = 5;
            else
                disp('fuck');
                cc(i) = 6;
            end
        end
        
        cell_classes = [cell_classes; cc];
        tuned_cells = [tuned_cells; is_tuned];
        s_classes{iFile} = cc(is_tuned);
    end
    
    if reassign_others
        mem_ind = abs(dpd_wo) ./ min( abs(dpd_ad) , abs(dpd_adwo) );
        idx = cell_classes == 5 & mem_ind >= 1;
        cell_classes(idx) = 3;
        idx = cell_classes == 5 & mem_ind < 1;
        cell_classes(idx) = 2;
    end
    
    tuned_cells = logical(tuned_cells);
    
    dpd_ad = dpd_ad(tuned_cells);
    dpd_wo = dpd_wo(tuned_cells);
    dpd_adwo = dpd_adwo(tuned_cells);
    cell_classes = cell_classes(tuned_cells);
    
    [k,x] = hist(cell_classes,1:5);
    bar(x+0.25*(iBlock-1),100*k/sum(k),0.25)
    set(gca,'Box','off','TickDir','out','FontSize',14,'XTick',1.25:5.25,'XTickLabel',{'Kin','Dyn','MemI','MemII','Other'},'XLim',[0.75 6]);
    ylabel('Percent of Cells','FontSize',14);
    title([use_model '-based neurons'],'FontSize',14);
end


%% plot sliding window thing
if 1
    %     root_dir = 'F:\trial_data_files\biomech_sim_results\';
    use_model = 'muscle';
    doMD = false;
    doAbs = true;
    
    ymin = 0;
    ymax = 180;
    
    figure;
    subplot1(1,2);
            filenames = { ...
%             'Chewie_CO_FF_2013-10-22', ...
%             'Chewie_CO_FF_2013-10-23', ...
%             'Chewie_CO_FF_2013-10-31', ...
%             'Chewie_CO_FF_2013-11-01', ...
%             'Chewie_CO_FF_2013-12-03', ...
%                     'Chewie_CO_FF_2013-12-04', ...
%             'Chewie_CO_FF_2015-07-01', ...
            'Chewie_CO_FF_2015-07-03', ...
%             'Chewie_CO_FF_2015-06-29', ...
%             'Chewie_CO_FF_2015-06-30', ...
%             'Chewie_CO_FF_2015-07-06', ...
%                     'Chewie_CO_FF_2015-07-07', ...
%                             'Chewie_CO_FF_2015-07-08', ...
            };
    
    for iFile = 1:length(filenames)
        iFile
        load(fullfile(dataDir,[filenames{iFile} '_tuning.mat']),'sw_data','tc_data');
        
        %tuned_cells = mean(tc_data(1).(use_model).rs,2) > 0.5;
        tuned_cells = all(prctile(tc_data(1).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
            all(prctile(tc_data(2).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
            all(prctile(tc_data(3).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
            angleDiff(tc_data(1).(use_model).cb{3}(:,1),tc_data(1).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
            angleDiff(tc_data(2).(use_model).cb{3}(:,1),tc_data(2).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
            angleDiff(tc_data(3).(use_model).cb{3}(:,1),tc_data(3).(use_model).cb{3}(:,2),true,false) < 40*pi/180;
        
        if iFile == 1
            day_dpd = zeros(length(filenames),length(sw_data));
            all_dpd = cell(1,length(sw_data));
            all_bl = cell(1,length(sw_data));
            all_ad = cell(1,length(sw_data));
            all_f = zeros(length(filenames),length(sw_data));
        end
        
        for iBlock = 1:length(sw_data)
            if doMD
                bl = sw_data(iBlock).tc_data(1).(use_model).tc(tuned_cells,2);
                ad = sw_data(iBlock).tc_data(2).(use_model).tc(tuned_cells,2);
                dpd = abs(ad - bl);
                day_dpd(iFile,iBlock) = mean(dpd);
            else
                bl = sw_data(iBlock).tc_data(1).(use_model).tc(tuned_cells,3);
                ad = sw_data(iBlock).tc_data(2).(use_model).tc(tuned_cells,3);
                dpd = angleDiff(bl,ad,true,~doAbs);
                day_dpd(iFile,iBlock) = circular_mean(dpd);
            end
            
            df = abs(mean(sw_data(iBlock).data(2).f) - mean(sw_data(iBlock).data(1).f));
            
            all_dpd{iBlock} = [all_dpd{iBlock}; dpd];
            all_bl{iBlock} = [all_bl{iBlock}; bl];
            all_ad{iBlock} = [all_bl{iBlock}; ad];
            all_f(iFile,iBlock) = mean(df);
        end
    end
    
    all_dpd = cell2mat(all_dpd);
    all_bl = cell2mat(all_bl);
    all_ad = cell2mat(all_ad);
    
    subplot1(1);
    ax1 = gca;
    hold all;
    
    [m,s] = deal(zeros(1,size(all_dpd,2)));
    for i = 1:size(all_dpd,2)
        if doMD
            m(i) = mean(all_dpd(:,i));
            s(i) = std(all_dpd(:,i))./sqrt(length(all_dpd));
        else
            m(i) = circular_mean((all_dpd(:,i)));
            s(i) = circular_std(all_dpd(:,i))./sqrt(length(all_dpd));
        end
    end
    
    plot(1:size(all_f,2),m.*(180/pi),'b','LineWidth',2);
    plot([1:size(all_f,2);1:size(all_f,2)],[m-s; m+s].*(180/pi),'b','LineWidth',2);
    set(ax1,'XLim',[0 length(all_f)+1],'Box','off','TickDir','out','FontSize',14,'XTickLabel',[],'YLim',[ymin ymax]);
    ylabel('mean dPD','FontSize',14);
    xlabel('Windows over movement','FontSize',14);
    
    ax1_pos = get(ax1,'Position'); % store position of first axes
    ax2 = axes('Position',ax1_pos,...
        'YAxisLocation','right',...
        'Color','none', ...
        'TickDir','out');
    hold all;
    
    m = mean(all_f,1);
    s = std(all_f,1)./sqrt(size(all_f,1));
    plot(1:size(all_f,2),m,'k','LineWidth',2);
    plot([1:size(all_f,2);1:size(all_f,2)],[m-s; m+s],'k','LineWidth',2);
    set(ax2,'XLim',[0 length(all_f)+1],'Box','off','XTick',[],'FontSize',14,'YLim',[0 2.5],'YTick',[]);
    % ylabel('mean dForce','FontSize',14);
    
    title('Chewie','FontSize',14);
    
%         figure;
%         plot(reshape(all_f,1,size(all_f,2)*size(all_f,1)),reshape(day_dpd,1,size(day_dpd,2)*size(day_dpd,1)).*(180/pi),'bo','LineWidth',2);
%         monkey_dpd = reshape(day_dpd,1,size(day_dpd,2)*size(day_dpd,1)).*(180/pi);
%         monkey_f = reshape(all_f,1,size(all_f,2)*size(all_f,1));
%         drawnow;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     for i = 1:2:15
    %         figure;
    %         hist(all_dpd(:,i),-pi:pi/32:pi);
    %         pause;
    %         close all;
    %     end
    %
    %     keyboard;
    filenames = { ...
        'Mihili_CO_FF_2014-02-17', ...
        'Mihili_CO_FF_2014-02-18', ...
        'Mihili_CO_FF_2014-03-07', ...
        'Mihili_CO_FF_2015-06-11', ...
        'Mihili_CO_FF_2014-02-03', ...
        'Mihili_CO_FF_2015-06-17', ...
        'Mihili_CO_FF_2015-06-10', ...
        'Mihili_CO_FF_2015-06-15', ...
        'Mihili_CO_FF_2015-06-16', ...
        };
    
    for iFile = 1:length(filenames)
        iFile
        load(fullfile(dataDir,[filenames{iFile} '_tuning.mat']),'sw_data','tc_data');
        
        tuned_cells = all(prctile(tc_data(1).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
            all(prctile(tc_data(2).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
            all(prctile(tc_data(3).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
            angleDiff(tc_data(1).(use_model).cb{3}(:,1),tc_data(1).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
            angleDiff(tc_data(2).(use_model).cb{3}(:,1),tc_data(2).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
            angleDiff(tc_data(3).(use_model).cb{3}(:,1),tc_data(3).(use_model).cb{3}(:,2),true,false) < 40*pi/180;
        
        if iFile == 1
            day_dpd = zeros(length(filenames),length(sw_data));
            all_dpd = cell(1,length(sw_data));
            all_f = zeros(length(filenames),length(sw_data));
        end
        
        for iBlock = 1:length(sw_data)
            if doMD
                bl = sw_data(iBlock).tc_data(1).(use_model).tc(tuned_cells,2);
                ad = sw_data(iBlock).tc_data(2).(use_model).tc(tuned_cells,2);
                dpd = abs(ad - bl);
                day_dpd(iFile,iBlock) = mean(dpd);
            else
                bl = sw_data(iBlock).tc_data(1).(use_model).tc(tuned_cells,3);
                ad = sw_data(iBlock).tc_data(2).(use_model).tc(tuned_cells,3);
                dpd = angleDiff(bl,ad,true,~doAbs);
                day_dpd(iFile,iBlock) = circular_mean(dpd);
            end
            
            df = abs(mean(sw_data(iBlock).data(2).f) - mean(sw_data(iBlock).data(1).f));
            
            all_dpd{iBlock} = [all_dpd{iBlock}; dpd];
            
            all_f(iFile,iBlock) = mean(df);
        end
    end
    
    all_dpd = cell2mat(all_dpd);
    
    subplot1(2);
    ax1 = gca;
    hold all;
    
    for i = 1:size(all_dpd,2)
        if doMD
            m(i) = mean(all_dpd(:,i));
            s(i) = std(all_dpd(:,i))./sqrt(length(all_dpd));
        else
            m(i) = circular_mean(all_dpd(:,i));
            s(i) = circular_std(all_dpd(:,i))./sqrt(length(all_dpd));
        end
    end
    
    plot(1:size(all_f,2),m.*(180/pi),'b','LineWidth',2);
    plot([1:size(all_f,2);1:size(all_f,2)],[m-s; m+s].*(180/pi),'b','LineWidth',2);
    set(ax1,'XLim',[0 length(all_f)+1],'Box','off','TickDir','out','FontSize',14,'XTickLabel',[],'YTick',[],'YLim',[ymin ymax]);
    xlabel('Windows over movement','FontSize',14);
    
    ax1_pos = get(ax1,'Position'); % store position of first axes
    ax2 = axes('Position',ax1_pos,...
        'YAxisLocation','right',...
        'Color','none', ...
        'TickDir','out');
    hold all;
    
    
    m = mean(all_f,1);
    s = std(all_f,1)./sqrt(size(all_f,1));
    plot(1:size(all_f,2),m,'k','LineWidth',2);
    plot([1:size(all_f,2);1:size(all_f,2)],[m-s; m+s],'k','LineWidth',2);
    set(ax2,'XLim',[0 length(all_f)+1],'Box','off','XTick',[],'FontSize',14,'YLim',[0 2.5]);
    ylabel('mean dForce','FontSize',14);
    
    title('Mihili','FontSize',14);
    
    
    % temp_dpd = reshape(day_dpd,1,size(day_dpd,2)*size(day_dpd,1)).*(180/pi);
    % temp_f = reshape(all_f,1,size(all_f,2)*size(all_f,1));
    %
    % temp_f = temp_f(temp_dpd > 0);
    % temp_dpd = temp_dpd(temp_dpd > 0);
    %
    % hold all;
    %     plot(temp_f,temp_dpd,'ro','LineWidth',2);
    %     set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 3],'YLim',[30 150]);
    %     xlabel('Change in Force','FontSize',14);
    %     ylabel('dPD','FontSize',14);
    %     monkey_dpd = [monkey_dpd, temp_dpd];
    %     monkey_f = [monkey_f, temp_f];
    %     % fit regression line
    %     [b,bint,~,~,stats] = regress(monkey_dpd',[ones(length(monkey_f),1) monkey_f']);
    %     stats
    %     plot(monkey_f,(b(1)+b(2)*monkey_f),'k','LineWidth',2);
    
end
