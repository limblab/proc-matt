%%
K = 0.1*0.15;
% TH_c = 0;

td = trial_data(getTDidx(trial_data,'epoch','AD'));
f = zeros(length(td),2);
for iTrial = 1:length(td)
    v = td(iTrial).vel / 100;
    a = td(iTrial).acc / 100;
    
    Fc = zeros(size(v,1),2);
    for t = 1:size(v,1)
        Fc(t,:) = 100 * K * [cos(TH_c) * v(t,1) + sin(TH_c)*v(t,2), -sin(TH_c)*v(t,1) + cos(TH_c) * v(t,2)];
    end
    
    f(iTrial,:) = [rms(hypot(Fc(:,1),Fc(:,2))), rms((M1+M2)*hypot(a(:,1),a(:,2)))];
end

figure;
hist(f(:,1)./f(:,2))

%%
sd = sim_data(getTDidx(trial_data,'epoch','AD'));

t = zeros(length(sd),2);
for iTrial = 1:length(sd)
    Tp = sum(abs(sd(iTrial).torques_plan),2);
    Tf = sum(abs(sd(iTrial).torques_force),2);
    t(iTrial,:) = [rms(Tf(:,1)), rms(Tp(:,1))];
end

figure;
hist(t(:,1)./t(:,2));


%%
b=regress(t(:,1),[ones(size(t,1),1) f(:,1)]);
figure;
subplot(121);hold all;
plot(f(:,1),t(:,1),'d')
plot(f(:,1),b(1)+f(:,1)*b(2),'k-','LineWidth',3)
title(['b = ' num2str(b(2))])

b=regress(t(:,2),[ones(size(t,1),1) f(:,2)]);
subplot(122);hold all;
plot(f(:,2),t(:,2),'d')
plot(f(:,2),b(1)+f(:,2)*b(2),'k-','LineWidth',3)
title(['b = ' num2str(b(2))])