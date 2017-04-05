clear; clc; close all;

monkey = 'Chewie';
usedate = '2013-10-22';
% load tuning data and find tuned cells
load(['F:\results\m1_cf_JNS_results\' monkey '\' usedate '\M1_tuning\CO_FF_movement_regression_onpeak_' usedate '.mat']);

% load BL data
bl = load(['F:\results\m1_cf_JNS_results\' monkey '\' usedate '\CO_FF_BL_' usedate '.mat'],'M1');
% load AD data
ad = load(['F:\results\m1_cf_JNS_results\' monkey '\' usedate '\CO_FF_AD_' usedate '.mat'],'M1');
% load WO data
wo = load(['F:\results\m1_cf_JNS_results\' monkey '\' usedate '\CO_FF_WO_' usedate '.mat'],'M1');

%%
epochs = [1,4,7];
num_wf = 200;
is_tuned = find(all(classes.istuned,2));% & classes.classes(:,1) == 2);

for i = is_tuned'
        figure;
        % now plot tuning curves
    for j = 1:length(epochs)
        pd = tuning(epochs(j)).pds(i,:);
        % get firing rate to each target
        u = unique(tuning(epochs(j)).theta)';
        fr = zeros(1,length(u));
        for k = 1:length(u)
            fr(k) = nanmean(tuning(epochs(j)).fr(tuning(epochs(j)).theta==u(k),i));
        end
        
        u = [u,u(1)];
        fr = [fr,fr(1)];
        
        subplot(2,3,3+j);
        polar(u,fr,'k-');
        hold on;
        polar([pd(1),pd(1)],[0 max(fr)],'r-');
        polar([pd(2),pd(2)],[0 max(fr)],'r--');
        polar([pd(3),pd(3)],[0 max(fr)],'r--');
    end
    title(classes.classes(i));
    

    % plot bl waveforms
    subplot(2,3,1);
    temp = randperm(length(bl.M1.units(i).ts));
    plot(bl.M1.units(i).wf(:,temp(1:num_wf)),'k');
    title(bl.M1.units(i).id);
    % plot ad waveforms
    subplot(2,3,2);
    temp = randperm(length(ad.M1.units(i).ts));
    plot(ad.M1.units(i).wf(:,temp(1:num_wf)),'k');
    
    % plot wo waveforms
    subplot(2,3,3);
    temp = randperm(length(wo.M1.units(i).ts));
    plot(wo.M1.units(i).wf(:,temp(1:num_wf)),'k');

    min_y = Inf;
    max_y = -Inf;
    for j = 1:3
        subplot(2,3,j);
        V = axis;
        min_y = min([min_y, V(3:4)]);
        max_y = max([max_y, V(3:4)]);
    end
    for j = 1:3
        subplot(2,3,j);
        set(gca,'Box','off','TickDir','out','FontSize',12,'YLim',[min_y,max_y]);
    end

    pause;
    close all;
end


% for i = is_tuned'
%     figure;
%     subplot(2,2,1);
%     temp = randperm(length(bl.M1.units(i).ts));
%     plot(bl.M1.units(i).wf(:,temp(1:100)),'k');
%     set(gca,'XLim',[1 48],'YLim',[-1000,1000]);
%     subplot(2,2,2);
%     temp = randperm(length(ad.M1.units(i).ts));
%     plot(ad.M1.units(i).wf(:,temp(1:100)),'k');
%     set(gca,'XLim',[1 48],'YLim',[-1000,1000]);
%     
%     pds = zeros(3,7);
%     for j = 1:7
%         pds(:,j) = tuning(j).pds(i,:)*180/pi;
%     end
%     subplot(2,2,3:4); hold all;
%     plot(pds(1,:),'k');
%     plot(pds(2,:),'k--');
%     plot(pds(3,:),'k--');
%     pause;
%     close all;
% end