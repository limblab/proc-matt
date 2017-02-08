% load('F:\TrialDataFiles\Mihili_CO_FF_2014-02-17.mat')
load('F:\TrialDataFiles\Mihili_CO_FF_2015-06-15.mat')

%%
t_max = 111;

% get baseline speed
td = trial_data(getTDidx(trial_data,'epoch','BL'));

bl_speed = zeros(length(td),t_max);
for trial = 1:length(td)
    s = hypot(td(trial).vel(:,1),td(trial).vel(:,2));
    idx = td(trial).idx_peak_speed-40:td(trial).idx_peak_speed+70;
    s = s(idx);
    bl_speed(trial,:) = s;%interp1(1:length(s),s,linspace(1,length(s),1000));
end


% get three AD speeds
td_idx = find(getTDidx(trial_data,'epoch','AD'));

td = trial_data(td_idx(1:floor(0.33*length(td_idx))));
ad_speed1 = zeros(length(td),t_max);
for trial = 1:length(td)
    s = hypot(td(trial).vel(:,1),td(trial).vel(:,2));
    idx = td(trial).idx_peak_speed-40:td(trial).idx_peak_speed+70;
    s = s(idx);
    ad_speed1(trial,:) = s;%interp1(1:length(s),s,linspace(1,length(s),1000));
end

% td = trial_data(td_idx(ceil(0.33*length(td_idx)):floor(0.66*length(td_idx))));
% ad_speed2 = zeros(length(td),1000);
% for trial = 1:length(td)
%     s = hypot(td(trial).vel(:,1),td(trial).vel(:,2));
%     idx = td(trial).idx_movement_on-10:td(trial).idx_trial_end-10;
%     s = s(idx);
%     ad_speed2(trial,:) = interp1(1:length(s),s,linspace(1,length(s),1000));
% end

td = trial_data(td_idx(ceil(0.66*length(td_idx)):end));
ad_speed3 = zeros(length(td),t_max);
for trial = 1:length(td)
    s = hypot(td(trial).vel(:,1),td(trial).vel(:,2));
    idx = td(trial).idx_peak_speed-40:td(trial).idx_peak_speed+70;
    s = s(idx);
    ad_speed3(trial,:) = s;%interp1(1:length(s),s,linspace(1,length(s),1000));
end

%% plot things
s_max = 30;


figure;
subplot1(1,3);

subplot1(1);
s = bl_speed;
m = mean(s,1);
s = [m-std(s,[],1)./sqrt(size(s,1)); m+std(s,[],1)./sqrt(size(s,1))];
patch([1:t_max,t_max:-1:1], [s(2,:),fliplr(s(1,:))],[0.7,0.7,0.7])
plot(1:t_max,m,'k','LineWidth',2);
set(gca,'Box','off','TickDir','out','FontSize',14,'Xlim',[1 t_max],'YLim',[0 s_max])
ylabel('Speed (cm/s','FontSize',14);
title('Baseline','FontSize',14);

subplot1(2);
s = ad_speed1;
m = mean(s,1);
s = [m-std(s,[],1)./sqrt(size(s,1)); m+std(s,[],1)./sqrt(size(s,1))];
patch([1:t_max,t_max:-1:1], [s(2,:),fliplr(s(1,:))],[0.7,0.7,0.7])
plot(1:t_max,m,'k','LineWidth',2);
set(gca,'Box','off','TickDir','out','FontSize',14,'Xlim',[1 t_max],'YLim',[0 s_max])
title('Early Adaptation','FontSize',14);

subplot1(3);
s = ad_speed3;
m = mean(s,1);
s = [m-std(s,[],1)./sqrt(size(s,1)); m+std(s,[],1)./sqrt(size(s,1))];
patch([1:t_max,t_max:-1:1], [s(2,:),fliplr(s(1,:))],[0.7,0.7,0.7])
plot(1:t_max,m,'k','LineWidth',2);
set(gca,'Box','off','TickDir','out','FontSize',14,'Xlim',[1 t_max],'YLim',[0 s_max])
title('Late Adaptation','FontSize',14);
