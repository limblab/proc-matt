%% setup and config
clear; clc; close all;

monkey = 'Han';
use_date = '2017-02-06';

rewrite_cds = false;
max_vel = 50;
% 'R',82: reward (success)
% 'A',65: abort (fails trial pre-go-cue, e.g. leaves center target during delay period)
% 'F',70: fail (fails trial post-go-cue, e.g. doesn't make it to the target in time)
% 'I',74: incomplete (fails trial after entering outer target, e.g. doesn't hold)
tdInputArgs = struct( ...
    'trialResults',{{'R','F','I'}}, ... % which to include
    'excludeUnits',[255], ... % sort codes to exclude
    'binSize',0.01, ... % binning size in s
    'extraTime',[0.1 0.01]); % time before targ pres and after end in s

meta.angle_dir = 'CCW';
meta.rotation_angle = -1;
meta.force_magnitude = 0.15;
meta.force_angle = 1.48;
meta.perturbation = 'FF';

% kind of a lookup table for raeed's filenames and my epoch codes
epoch_names = {'baseline','BL'; 'baseline','BL';'adaptation','AD'; 'washout','WO'};

%% process to CDS and make trial data

trial_data = [];
for e = 1:4
    fin = fullfile('F:\',monkey,'CerebusData',use_date,[monkey '_' datestr(use_date,'yyyymmdd') '_CObumpcurl_' epoch_names{e,1} '_area2EMG_00' num2str(e) '.nev']);
    fout = fullfile('F:\',monkey,'CDS',use_date,[monkey '_' datestr(use_date,'yyyymmdd') '_CObumpcurl_' epoch_names{e,1} '_area2EMG_00' num2str(e) '.mat']);
    if ~exist(fout,'file') || rewrite_cds
%         cds=commonDataStructure();
%         cds.file2cds(fin,6,'arrayS1','monkeyHan','taskCObump','ranByRaeed','mapFile\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\L_MillerLab\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp','ignoreJumps')
%         % fname = 'F:\Chewie\CerebusData\2016-10-21\Chewie_PMd_RT_CS_BL_10212016_002.nev';
%         % cds.file2cds(fname,6,'arrayPMd','monkeyChewie','taskRW','ranByMatt','mapFileZ:\lab_folder\Animal-Miscellany\Chewie 8I2\Chewie Left PMd SN 6251-001469.cmp','ignoreJumps')
%         save(fout,'cds','-v7.3')
    else
        % load the existing CDS
        disp('CDS found, loading...');
        load(fout);
    end
    meta.epoch = epoch_names{e,2};
    tdInputArgs.meta = meta;
    td = parseFileByTrial_cds(cds,tdInputArgs);
    trial_data = [trial_data, td];
    clear cds;
end

save(fullfile('F:\',monkey,[monkey '_' datestr(use_date,'yyyymmdd') '_CObumpcurl.mat']),'trial_data');



%%
close all;

learning_params = struct('time_window',{{'idx_go_cue',20;'idx_trial_end',0}}, ...
    'result_codes',{{'R'}}, ...
    'use_bl_ref',false);
td = trial_data(getTDidx(trial_data,'result','R'));
m = get_learning_metrics(td,'corr',learning_params);

% learning_params = struct('time_window',{{'idx_go_cue',20;'idx_go_cue',40}}, ...
%     'result_codes',{{'R'}}, ...
%     'use_bl_ref',false);
% td = trial_data(getTDidx(trial_data,'epoch','BL','result','R'));
% m1 = get_learning_metrics(td(1:150),'angle',learning_params);
% 
% td = trial_data(getTDidx(trial_data,'epoch','AD','result','R'));
% m2 = get_learning_metrics(td(1:250),'angle',learning_params);
% 
% td = trial_data(getTDidx(trial_data,'epoch','WO','result','R'));
% m3 = get_learning_metrics(td(1:150),'angle',learning_params);
% m = [m1;m2;m3];

bins = 1:1:length(m)+1;
m_avg = zeros(1,length(bins)-1);
for i = 1:length(bins)-1
    m_avg(i) = mean(m(bins(i):bins(i+1)-1));
end
% 
% win_size = 3;
% m_avg = zeros(1,length(m)-win_size);
% for i = 1:length(m)-win_size
%     m_avg(i) = mean(m(i:i+win_size));
% end

figure; hold all;
% plot([15,15],[-1.5,1.5]*180/pi,'k--','LineWidth',2);
% plot([30,30],[-1.5,1.5]*180/pi,'k--','LineWidth',2);
plot(m_avg,'+');
axis('tight');
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[-1.5,1.5]);
