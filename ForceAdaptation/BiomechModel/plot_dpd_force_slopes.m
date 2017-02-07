clear;
close all;
clc;

dataSummary;
dataDir = fullfile(rootDir,resultsDir,'biomech_model');

use_model = 'muscle';
doMD = false;
doAbs = true;

min_r2 = 0.5;
ymin = 0;
ymax = 100;

% filenames = { 'tuning_01.mat','tuning_05.mat','tuning_025.mat','tuning_1.mat','tuning_2.mat','tuning_10.mat'};
% filenames = {'tuning_120.mat','tuning_12.mat','tuning_1.mat'};


filenames = {'tuning_001.mat','tuning_005.mat','tuning_01.mat','tuning_015.mat','tuning_02.mat','tuning_025.mat','tuning_05','tuning_075','tuning_1.mat','tuning_2.mat'};
plot_colors = distinguishable_colors(length(filenames));


figure;
for iFile = 1:length(filenames)
    load(fullfile(dataDir,'changing_beta',filenames{iFile}),'sw_data','tc_data');
    
    %tuned_cells = mean(tc_data(1).(use_model).rs,2) > 0.5;
    tuned_cells = all(prctile(tc_data(1).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
        all(prctile(tc_data(2).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
        all(prctile(tc_data(3).(use_model).rs,[2.5 97.5],2) > min_r2,2) & ...
        angleDiff(tc_data(1).(use_model).cb{3}(:,1),tc_data(1).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
        angleDiff(tc_data(2).(use_model).cb{3}(:,1),tc_data(2).(use_model).cb{3}(:,2),true,false) < 40*pi/180 & ...
        angleDiff(tc_data(3).(use_model).cb{3}(:,1),tc_data(3).(use_model).cb{3}(:,2),true,false) < 40*pi/180;

        day_dpd = zeros(length(filenames),length(sw_data));
        all_dpd = cell(1,length(sw_data));
        all_bl = cell(1,length(sw_data));
        all_ad = cell(1,length(sw_data));
        all_f = zeros(length(filenames),length(sw_data));
    
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
    
    all_dpd = cell2mat(all_dpd);
    all_bl = cell2mat(all_bl);
    all_ad = cell2mat(all_ad);
    
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
    
    subplot(131);
    hold all;
    plot(1:size(all_f,2),m.*(180/pi),'LineWidth',2,'Color',plot_colors(iFile,:));
    plot([1:size(all_f,2);1:size(all_f,2)],[m-s; m+s].*(180/pi),'LineWidth',2,'Color',plot_colors(iFile,:));
    
    subplot(132);
    hold all;
    plot(iFile,180/pi*(m(end)-m(1))/5,'o','LineWidth',3,'Color',plot_colors(iFile,:));
    
end

subplot(131);
set(gca,'XLim',[0 7],'Box','off','TickDir','out','FontSize',14,'XTickLabel',[]);
ylabel('mean dPD','FontSize',14);
xlabel('Windows over movement','FontSize',14);

subplot(132);
set(gca,'XLim',[0 length(filenames)+1],'Box','off','TickDir','out','FontSize',14,'XTickLabel',[]);
ylabel('Slope','FontSize',14);
xlabel('Windows over movement','FontSize',14);

subplot(133);
hold all;
for iFile = 1:length(filenames)
    plot(0,0,'o','Color',plot_colors(iFile,:),'LineWidth',5);
end
set(gca,'Box','off','TickDir','out','FontSize',14,'Visible','off');
legend({'0.01','0.05','0.1','0.15','0.2','0.25','0.5','0.75','1','2'},'FontSize',14);
