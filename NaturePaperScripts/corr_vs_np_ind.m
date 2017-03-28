%% Prep datasets and get paths etc
clear; clc; close all;
dataSummary;

% Session parameters

task = 'CO';
pert = 'FF';

sessions = { ...
    'Chewie','2016-09-15'; ... % CF
    'Chewie','2016-09-21'; ...
    'Chewie','2016-10-05'; ...
    'Chewie','2016-10-07'; ...
    'Chewie','2016-10-11'; ...
    'Mihili','2014-02-03'; ...
    'Mihili','2014-02-17'; ...
    'Mihili','2014-02-18'; ...
    'Mihili','2014-03-07'; ...
    };
monkeys = unique(sessions(:,1));
dates = unique(sessions(:,2));

session_idx = getFileDBidx(filedb, ...
    {'Task',task,'Perturbation',pert,'Monkey',monkeys,'Date',dates});

idx_start = {'idx_target_on', 0};
idx_end   = {'idx_trial_end', -20};

badtrial_params = struct(...
    'ranges', {{'idx_go_cue','idx_movement_on',[5 40]; ...
    'idx_peak_speed','idx_trial_end',[60 100]}});

badneuron_params_m1 = struct( ...
    'arrays','M1', ...
    'min_fr',1, ...
    'do_shunt_check',0, ...
    'use_trials',{{'epoch',{'BL','AD'}}});

badneuron_params_pmd = struct( ...
    'arrays','PMd', ...
    'min_fr',1, ...
    'do_shunt_check',0, ...
    'use_trials',{{'epoch',{'BL','AD'}}});

func_calls = { ...
    {@getTDidx,'result','R'}, ...
    {@removeBadTrials,badtrial_params}, ...
    {@removeBadNeurons,badneuron_params_m1}, ...
    {@removeBadNeurons,badneuron_params_pmd}, ...
    @sqrtTransform, ...
    {@smoothSignals,struct('signals',{{'M1_spikes','PMd_spikes'}},'calc_fr',true,'kernel_SD',0.1)}, ...
    {@trimTD,idx_start,idx_end}};


[all_dr_m1,all_dr_pmd,all_dr_m1pmd,all_np_ind] = deal([]);
% load it
figure; hold all;
for file = session_idx'
    fname = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
    [td,params] = loadTDfiles(fname, func_calls{:});
    
    [~,td1] = getTDidx(td,'epoch','BL');
    [~,td2] = getTDidx(td,'epoch','AD','range',[0.33 1]);
    td_s = [td1, td2];
    
    %td = appendTDs(trimTD(td_s,{'idx_target_on',0},{'idx_target_on',50}),trimTD(td_s,{'idx_go_cue',-0},{'idx_go_cue',40}));
    td = trimTD(td_s,{'idx_go_cue',-10},{'idx_go_cue',60});
    td = binTD(td,3);
    % td = softNormalize(td);
    td = trialAverage(td,{'target_direction','epoch'});
    
    % [~,td0] = getTDidx(td,'epoch','BL','range',[0 0.66]);
    [~,td1] = getTDidx(td,'epoch','BL','range',[0 1]);
    [~,td2] = getTDidx(td,'epoch','AD','range',[0.33 1]);
    
     
    arrays = {'M1_spikes'};
    params_corr = struct('signals',{arrays},'cluster_order',false,'do_norm',false,'cluster_arrays',false);
    [m1_bl,m1_bl_sort] = pairwiseCorr(td1,params_corr);
    [m1_ad,m1_ad_sort] = pairwiseCorr(td2,params_corr);
    arrays = {'PMd_spikes'};
    params_corr = struct('signals',{arrays},'cluster_order',false,'do_norm',false,'cluster_arrays',false);
    [pmd_bl,pmd_bl_sort] = pairwiseCorr(td1,params_corr);
    [pmd_ad,pmd_ad_sort] = pairwiseCorr(td2,params_corr);
    arrays = {'M1_spikes','PMd_spikes'};
    params_corr = struct('signals',{arrays},'cluster_order',false,'do_norm',false,'cluster_arrays',true);
    [m1pmd_bl,m1pmd_bl_sort] = pairwiseCorr(td1,params_corr);
    [m1pmd_ad,m1pmd_ad_sort] = pairwiseCorr(td2,params_corr);
    m1pmd_bl = m1pmd_bl(1:size(m1_bl,1),size(m1_bl,1)+1:end);
    m1pmd_ad = m1pmd_ad(1:size(m1_bl,1),size(m1_bl,1)+1:end);
    
    
    % get null and potent spaces
    %td = appendTDs(trimTD(td,{'idx_target_on',0},{'idx_target_on',50}),trimTD(td_s,{'idx_go_cue',-0},{'idx_go_cue',40}));
    td = appendTDs( ...
            trimTD(td_s,{'idx_target_on',0},{'idx_target_on',50}), ...
            trimTD(td_s,{'idx_go_cue',0},{'idx_go_cue',20}), ...
            trimTD(td_s,{'idx_movement_on',-10},{'idx_movement_on',50}));
    
    pn_params = struct( ...
        'in_signals','PMd_spikes', ...
        'out_signals','M1_spikes', ...
        'in_dims',16, ...
        'out_dims',8, ...
        'use_trials',getTDidx(td,'epoch','BL'));
    [td,pn_info] = getPotentSpace(td,pn_params);
    
    % get weight of each PMd neuron onto the null and potent spaces
    w_potent = pn_info.w_in(:,1:16)*pn_info.V_potent;
    w_null = pn_info.w_in(:,1:16)*pn_info.V_null;
    
    % get change in correlation after larning
    dr = pmd_ad - pmd_bl;
    for i = 1:size(dr,1), dr(i,i) = NaN; end
    dr = nanmean(dr,2);
    
    % define null/potent index
    np_ind = (rms(w_potent,2)-rms(w_null,2))./(rms(w_potent,2)+rms(w_null,2));
    
    all_np_ind = [all_np_ind; np_ind];
    all_dr_pmd = [all_dr_pmd; dr];
    
    dr = m1_ad - m1_bl;
    for i = 1:size(dr,1), dr(i,i) = NaN; end
    dr = nanmean(dr,2);
    all_dr_m1 = [all_dr_m1; dr];
    
    dr = m1pmd_ad - m1pmd_bl;
    for i = 1:size(dr,1), dr(i,i) = NaN; end
    dr = nanmean(dr,2);
    all_dr_m1pmd = [all_dr_m1pmd; dr];
end

plot(all_np_ind, all_dr_pmd, '.','LineWidth',2);
axis('square');
set(gca,'Box','off','TickDir','out','FontSize',14);

V = axis;
x_min = -max(abs(V(1:2)));
x_max = max(abs(V(1:2)));
y_min = -max(abs(V(3:4)));
y_max = max(abs(V(3:4)));

plot([0 0],[y_min,y_max],'k--');
plot([x_min,x_max],[0 0]);
set(gca,'XLim',[x_min x_max],'YLim',[y_min y_max]);
xlabel('Null/Potent Index','FontSize',14);
ylabel('Change in correlation','FontSize',14);

%%
figure;
hist(all_np_ind,-0.6:0.05:0.6);
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[-0.6 0.6]);
ylabel('Count');
xlabel('<--- MORE NULL     Potent/Null Index     MORE POTENT--->');


