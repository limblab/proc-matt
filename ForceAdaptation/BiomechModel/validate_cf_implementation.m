%%
TH_c = 85*pi/180;
K = 0.15;

v = -50+100*rand(100000,2);

Fc = zeros(size(v));
for t = 1:size(v,1)
    %Fc(t,:) = 100 * K * [cos(TH_c)*v(t,1) + sin(TH_c)*v(t,2), sin(TH_c)*v(t,1) - cos(TH_c) * v(t,2)];
    Fc(t,:) = 100 * K * [cos(TH_c)*v(t,1) + sin(TH_c)*v(t,2), -sin(TH_c)*v(t,1) + cos(TH_c) * v(t,2)];
end

close all;
figure;
hist(angleDiff(atan2(Fc(:,2),Fc(:,1)),atan2(v(:,2),v(:,1)),true,true).*180/pi)
figure;
plot(hypot(v(:,1),v(:,2)),hypot(Fc(:,1),Fc(:,2)),'.')

%%
TH_c = pi-85*pi/180;
K = 0.15;

close all;
figure;
hold all;
for v1 = -30:5:30
    for v2 = -30:5:30
        Fc = K * [cos(TH_c)*v1 + sin(TH_c)*v2, -sin(TH_c)*v1 + cos(TH_c) * v2];
%         arrow([v1,v2],[v1+Fc(1),v2+Fc(2)]);
    plot(v1,v2,'ko');
    plot([v1,v1+Fc(1)],[v2,v2+Fc(2)],'k');
    end
end

axis('tight');
axis('square');
set(gca,'Box','off','TickDir','out','FontSize',14,'XTick',[-30,-15,0,15,30],'YTick',[-30,-15,0,15,30]);
xlabel('X-Velocity','FontSize',14);
ylabel('Y-Velocity','FontSize',14);