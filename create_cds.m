%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load stuff
dataSummary;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to do
remakeCDS = true; % if false, will skip CDS generating if file exists
remakeTD = false; % make trial data format that groups epochs
remakeFileDB = false; % add CDS info to fileDB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filedb = filedb_add(filedb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% which session
% which_sessions = getFileDBidx(filedb,{'Date','2016-09-21'});
%     {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
%     'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});
which_sessions = find(strcmpi(filedb.Date,'2016-10-21') & strcmpi(filedb.Task,'CO'));
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
%     .extra_time    : [time before, time after] beginning and end of trial (default [0.2 0.2] sec)
%     .all_points    : flag to include all data points. Thus, disregards extra_time
%                       and each trial ends at trial_start of the one after
%     .pos_offset    : offset (in units of cds.pos) to zero position (default [0,0])
max_vel = 100;
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
    'extra_time',[0.1 0.1], ... % time before targ pres and after end in s
    'all_points',false);
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
    
    if isempty(s_idx)% if it's not in the database yet, add it
        filedb = filedb_add(filedb);
        s_idx = find(strcmpi(filedb.Date,date) & ...
            strcmpi(filedb.Monkey,monkey) & ...
            strcmpi(filedb.Task,task),1);
    end
    
    pert = filedb.Perturbation{s_idx};
    % get list of epochs and filenums
    % Note: assumes Ricardo's merge and split function was used for sorting
    % Find all of the .mat files
    d = dir(fullfile(rootDir,monkey,cerebusDataDir,date,'*.mat'));
    d = {d.name};
    d = d(cellfun(@(x) isempty(regexp(x,'-s.mat', 'once')) & ...
        isempty(regexp(x,'-metatags.mat', 'once')) & ...
        isempty(regexp(x,'_PMd_', 'once')) & ...
        ~isempty(regexp(x,['_' task '_'], 'once')),d));
    
    if ~isempty(d)
        % name format will be MONKEY_ARRAY_TASK_PERT_EPOCH_DATE_NUM.mat
        [filenums,I] = sort(cellfun(@(x) x(end-6:end-4),d,'uni',0));
        %   date is always 8 characters and num is always 3, so epoch is end-18:end-17
        epochs = cellfun(@(x) x(end-18:end-17),d,'uni',0);
        epochs = epochs(I); % resort to match filenums
        
        if str2double(date(1:4)) < 2016 % all Lab 3 data is 2015 or earlier
            lab = 3;
        else
            lab = 6;
        end
        
        % get array maps
        arrays = array_list{strcmpi(array_list(:,1),monkey) & cellfun(@(x) datenum(date) < datenum(x), array_list(:,2)),3};
        arrayMaps = array_list{strcmpi(array_list(:,1),monkey) & cellfun(@(x) datenum(date) < datenum(x), array_list(:,2)),4};
        
        td_filename = fullfile(rootDir,TDDir,[monkey '_' task '_' pert '_' date '.mat']);
        
        if ~exist(td_filename,'file') || remakeTD
            trial_data = [];
            meta.perturbation = pert;
            switch lower(filedb.Direction{s_idx})
                case 'ccw'
                    perturbation_direction = 1;
                case 'cw'
                    perturbation_direction = -1;
            end
            % for VR: angle (negative is clockwise); for FF: [magnitude, direction] (negative is clockwise)
            switch lower(meta.perturbation)
                case 'ff'
                    meta.perturbation_info = [filedb.FFMagnitude(s_idx), perturbation_direction*filedb.FFAngle(s_idx)];
                case 'vr'
                    meta.perturbation_info = perturbation_direction*filedb.VRAngle(s_idx);
            end
            
        end
        
        somethingCDS = false; %flag variable
        
        % loop along all of the files
        for f = 1:length(filenums)
            epoch = epochs{f};
            filenum = filenums{f};
            
            outname = fullfile(rootDir,monkey,CDSDir,date,[monkey '_' task '_' pert '_' epoch '_' datestr(date,'mmddyyyy') '_' filenum '.mat']);
            if ~exist(outname,'file') || remakeCDS
                somethingCDS = true; % we are doing something CDS related
                
                fname = fullfile(rootDir,monkey,cerebusDataDir,date,[monkey '_' arrays{1} '_' task '_' pert '_' epoch '_' datestr(date,'mmddyyyy') '_' filenum '.nev']);
                
                % make blank cds class:
                cds=commonDataStructure();
                % some older data had the wrong sync label for cerebus 1
                cds.aliasList = {'Cerebus2Sync','Cerebus1Sync'};
                %load the data:
                if strcmpi(task,'rt')
                    cds_task = 'RW';
                else
                    cds_task = task;
                end
                cds.file2cds(fname,lab,['array' arrays{1}],['monkey' monkey],['task' cds_task],'ranByMatt',['mapFile' arrayMaps{1}],'ignoreJumps')
                
                if length(arrays) > 1 % add in the PMd data
                    fname = fullfile(rootDir,monkey,cerebusDataDir,date,[monkey '_' arrays{2} '_' task '_' pert '_' epoch '_' datestr(date,'mmddyyyy') '_' filenum '.nev']);
                    cds.aliasList = {};
                    % load spiking for PMd array
                    cds.file2cds(fname,lab,['array' arrays{2}],['monkey' monkey],['task' cds_task],'ranByMatt',['mapFile' arrayMaps{2}],'ignoreJumps')
                end
                
                % save it
                if ~exist(fullfile(rootDir,monkey,CDSDir,date),'dir')
                    mkdir(fullfile(rootDir,monkey,CDSDir,date));
                end
                save(outname,'-v7.3','cds')
            elseif ~exist(td_filename,'file') || remakeTD
                disp([epoch ' CDS already exists... loading...']);
                load(outname,'cds');
            end
            
            if ~exist(td_filename,'file') || remakeTD
                disp('Converting to TrialData struct...');
                meta.epoch = epoch;
                tdInputArgs.meta = meta;
                td = parseFileByTrial(cds,tdInputArgs);
                
                trial_data = [trial_data, td];
            end
            clear td;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check quality of trials and prune some bad ones
        if ~exist(td_filename,'file') || remakeTD            
            disp('Checking sorted units...');
            trial_data = removeBadTrials(trial_data);
            trial_data = getCommonUnits(trial_data);
            
            if strcmpi(task,'rt')
                for i = 1:length(trial_data)
                    trial_data(i).task = 'RT';
                end
            else % can't do onset for RT yet
                trial_data = getMoveOnsetAndPeak(trial_data);
            end
            save(td_filename,'trial_data');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add CDS info to filedb
        if somethingCDS || remakeFileDB % assume if this was false that there is nothing new to add
            disp('Adding CDS data to FileDB...');
            filedb = filedb_add(filedb,fullfile(rootDir,monkey,CDSDir,date));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        warning(['File for ' date ' has not been sorted or split yet, apparently...']);
    end
end
