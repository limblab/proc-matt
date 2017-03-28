load('/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-07.mat')
% load('/Users/mattperich/Data/TrialDataFiles/Mihili_CO_FF_2014-02-03.mat')
[~,td] = getTDidx(trial_data,'result','R');

td = trimTD(td,{'idx_go_cue',0},{'idx_go_cue',60});
td = sqrtTransform(td);
td = smoothSignals(td,struct('signals',{{'M1_spikes','PMd_spikes'}},'kernel_SD',0.1));
pn_params = struct( ...
    'in_signals','PMd_spikes', ...
    'out_signals','M1_spikes', ...
    'in_dims',16, ...
    'out_dims',8, ...
    'use_trials',{{'epoch','BL'}});
td = getPotentSpace(td,pn_params);

[~,td_temp] = getTDidx(td,'epoch','BL');
[A,B,R_bl,U,V] = canoncorr(cat(1,td_temp.PMdM1_potent),cat(1,td_temp.PMdM1_null));
figure;
plot(R_bl,'LineWidth',3);
set(gca,'Box','off','TickDir','out','FontSize',14);


%%
[~,td_bl1] = getTDidx(td,'epoch','BL','range',[0 0.5]);
[~,td_bl2] = getTDidx(td,'epoch','BL','range',[0.5 1]);
[~,~,R_bl_p,~,~] = canoncorr(cat(1,td_bl1.PMdM1_potent),cat(1,td_bl2.PMdM1_potent));
[~,~,R_bl_n,~,~] = canoncorr(cat(1,td_bl1.PMdM1_null),cat(1,td_bl2.PMdM1_null));

[~,td_temp1] = getTDidx(td,'epoch','BL','range',[10 110]);
[~,td_temp2] = getTDidx(td,'epoch','AD','range',[100 200]);
[A,B,R_ad_p,U,V] = canoncorr(cat(1,td_temp1.PMdM1_potent),cat(1,td_temp2.PMdM1_potent));
[A,B,R_ad_n,U,V] = canoncorr(cat(1,td_temp1.PMdM1_null),cat(1,td_temp2.PMdM1_null));

[~,td_temp2] = getTDidx(td,'epoch','AD','range',[1 101]);
[A,B,R_ad_p2,U,V] = canoncorr(cat(1,td_temp1.PMdM1_potent),cat(1,td_temp2.PMdM1_potent));
[A,B,R_ad_n2,U,V] = canoncorr(cat(1,td_temp1.PMdM1_null),cat(1,td_temp2.PMdM1_null));

figure;
subplot(121); hold all;
plot(R_bl_p,'LineWidth',3);
plot(R_ad_p,'LineWidth',3);
plot(R_ad_p2,'LineWidth',3);
set(gca,'Box','off','TickDir','out','FontSize',14);
title('Potent Space');
legend({'Baseline','Early Learning','Late Learning'});
ylabel('Canonical Correlation')
xlabel('Dimensions')

subplot(122); hold all;
plot(R_bl_n,'LineWidth',3);
plot(R_ad_n,'LineWidth',3);
plot(R_ad_n2,'LineWidth',3);
set(gca,'Box','off','TickDir','out','FontSize',14);
title('Null Space');




