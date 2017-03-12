%%
targ_angs = [pi/2, pi/4, 0, -pi/4, -pi/2, -3*pi/4, pi, 3*pi/4];
a1 = atan2(cds.trials.tgtCtr(:,2),cds.trials.tgtCtr(:,1));
a2 = cds.trials.tgtDir;
% a3 = a2; a3(~isnan(a3)) = targ_angs(cds.trials.tgtID(~isnan(a3))+1);


%%
t = cds.kin.t;
for trial = find(cds.trials.tgtID == 12)%1:size(cds.trials,1)
    if strcmpi(cds.trials.result(trial),'r')
    t_start = cds.trials.tgtOnTime(trial);
    t_end = cds.trials.endTime(trial);
    
    idx = find(t >= t_start & t <= t_end);
    figure; hold all;
    plot(cds.kin.x(idx)-3,cds.kin.y(idx)+32,'k-');
    plot(cds.kin.x(idx(end))-3,cds.kin.y(idx(end))+32,'ko','LineWidth',3);
    plot(8*cos(a1(trial)),8*sin(a1(trial)),'bo','LineWidth',3);
    plot(8*cos(a2(trial)),8*sin(a2(trial)),'ro','LineWidth',3);
%     plot(8*cos(a3(trial)),8*sin(a3(trial)),'go','LineWidth',3);
    set(gca,'YLim',[-9,9],'XLim',[-9,9]);
    pause;
    close all;
    end
end