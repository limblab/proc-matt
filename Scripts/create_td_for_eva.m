%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load stuff
dataSummary;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
monkey = 'Chewie';
task = 'CO';
pert = 'VR';
date = '2016-10-06';
epoch = 'BL';
filenum = '001';
SubDir = 'eva';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
which_sessions = ismember(filedb.Task,task) & ...
    strcmpi(filedb.Monkey,monkey) & ...
    ismember(filedb.Perturbation,pert) & ...
    datenum(filedb.Date) == datenum(date);
% Chewie: 2013-12-04, 2013-12-19, 2015-07-07, 2015-07-09 had CDS bugs
% Chewie '2015-07-08' didn't have a washout in trial_table
sessions = [filedb.Monkey(which_sessions), filedb.Date(which_sessions), filedb.Task(which_sessions)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters for the trial_data
%     .meta          : a struct with a field for each meta parameter you want attached
%                       to this file. This can handle any arbitrary information!
%     .event_names   : Which cds.trials events to add to struct
%                       Format: {'CDS_TRIAL_TABLE_NAME','ALIAS'; ... etc ...}
%                       Can ignore ALIAS and give Nx1 cell vector if you want
%                       By default, assumes startTime and endTime exist,
%                       and attempts to add tgtOnTime and goCueTime if possible
%     .array_alias   : Aliases for renaming arrays from CDS names
%                       Format: {'CDS_NAME','NEW_NAME'; ...etc...}
%     .exclude_units : ID for which units to exclude (Default: [0,255])
%                       NOTE: this default gets rid of unsorted!
%     .trial_results : which reward codes to use ('R','A','F','I') (Default: {'R'})
%     .bin_size      : default 0.01 sec
%     .extra_time    : [time before, time after] beginning and end of trial (default [0.2 0.2] sec)
%     .all_points    : flag to include all data points. Thus, disregards extra_time
%                       and each trial ends at trial_start of the one after
%     .pos_offset    : offset (in units of cds.pos) to zero position (default [0,0])
trial_results = {'R','F','I'};
event_list = { ...
    'startTime','trial_start'; ...
    'tgtOnTime','target_on'; ...
    'goCueTime','go_cue'; ...
    'endTime','trial_end'};
% 'R',82: reward (success)
% 'A',65: abort (fails trial pre-go-cue, e.g. leaves center target during delay period)
% 'F',70: fail (fails trial post-go-cue, e.g. doesn't make it to the target in time)
% 'I',74: incomplete (fails trial after entering outer target, e.g. doesn't hold)
tdInputArgs = struct( ...
    'event_list',{event_list}, ...
    'trial_results',{trial_results}, ... % which to include
    'exclude_units',[0,255], ... % sort codes to exclude
    'bin_size',0.01, ... % binning size in s
    'extra_time',[0.1 0.1], ... % time before targ pres and after end in s
    'all_points',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build CDS and trial data
% if multiple dates are provided, loop along them
for iS = 1:size(sessions,1)
    disp(['Session ' num2str(iS) ' of ' num2str(size(sessions,1))]);
    monkey = sessions{iS,1};
    date = sessions{iS,2};
    task = sessions{iS,3};
    
    s_idx = find(strcmpi(filedb.Date,date) & strcmpi(filedb.Monkey,monkey) & strcmpi(filedb.Task,task),1);
    
    td_filename = fullfile(rootDir,TDDir,SubDir,[monkey '_' task '_' pert '_' date '.mat']);
    
    outname = fullfile(rootDir,monkey,CDSDir,date,[monkey '_' task '_' pert '_' epoch '_' datestr(date,'mmddyyyy') '_' filenum '.mat']);
    
    if ~exist(td_filename,'file') || remakeTD
        disp([epoch ' CDS already exists... loading...']);
        load(outname,'cds');
    end
    
    disp('Converting to TrialData struct...');
    trial_data = parseFileByTrial(cds,tdInputArgs);
end
clear td;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check quality of trials and prune some bad ones
disp('Checking sorted units...');
trial_data = getCommonUnits(trial_data);
save(td_filename,'trial_data');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

