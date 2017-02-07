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
    'excludeUnits',[0,255], ... % sort codes to exclude
    'binSize',0.01, ... % binning size in s
    'extraTime',[0.5 0.5]); % time before targ pres and after end in s

meta.angle_dir = 'CCW';
meta.rotation_angle = -1;
meta.force_magnitude = 0.15;
meta.force_angle = 1.48;
meta.perturbation = 'FF';

% kind of a lookup table for raeed's filenames and my epoch codes
epoch_names = {'baseline','BL'; 'adaptation','AD'; 'washout','WO'};

%% process to CDS and make trial data

trial_data = [];
for e = 1:3
    fin = fullfile('F:\',monkey,'CerebusData',use_date,[monkey '_' datestr('2017-01-26','yyyymmdd') '_CObumpcurl_' epoch_names{e,1} '_area2EMG_00' num2str(e) '.nev']);
    fout = fullfile('F:\',monkey,'CDS',use_date,[monkey '_' datestr('2017-01-26','yyyymmdd') '_CObumpcurl_' epoch_names{e,1} '_area2EMG_00' num2str(e) '.mat']);
    if ~exist(fout,'file') || rewrite_cds
        cds=commonDataStructure();
        cds.file2cds(fin,6,'arrayS1','monkeyHan','taskCObump','ranByRaeed','mapFile\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\L_MillerLab\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp','ignoreJumps')
        % fname = 'F:\Chewie\CerebusData\2016-10-21\Chewie_PMd_RT_CS_BL_10212016_002.nev';
        % cds.file2cds(fname,6,'arrayPMd','monkeyChewie','taskRW','ranByMatt','mapFileZ:\lab_folder\Animal-Miscellany\Chewie 8I2\Chewie Left PMd SN 6251-001469.cmp','ignoreJumps')
        save(fout,'cds','-v7.3')
    else
        % load the existing CDS
        disp('CDS found, loading...');
        load(fout);
    end
    meta.epoch = epoch_names{e,2};
    tdInputArgs.meta = meta;
    td = parseFileByTrial_cds(cds,tdInputArgs);
    trial_data = [trial_data, td];
end

save(fullfile('F:\',monkey,[monkey '_' datestr('2017-01-26','yyyymmdd') '_CObumpcurl.mat']),'trial_data');

%%
close all;
td = trial_data(~isnan([trial_data.target_direction]));
m = get_learning_metrics(td(getTDidx(td,'epoch','AD')),'angle');

m(m < -0.2) = [];
plot(m,'o');
