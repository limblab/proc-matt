%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load stuff
dataSummary;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to do
remakeCDS = false; % if false, will skip CDS generating if file exists
makeTD = false; % make trial data format that groups epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filedb = filedb_add(filedb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% which session
which_sessions = strcmpi(filedb.Monkey,'Chewie') & ...
                 strcmpi(filedb.Task,'RT') & ...
                 ismember(filedb.Perturbation,{'CS'}) & ...
                 datenum(filedb.Date) == datenum('2016-10-21');
             % Chewie: 2013-12-04, 2013-12-19, 2015-07-07, 2015-07-09 had CDS bugs
             % Chewie '2015-07-08' didn't have a washout in trial_table
sessions = [filedb.Monkey(which_sessions), filedb.Date(which_sessions)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters for the trial_data
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build CDS and trial data

% if multiple dates are provided, loop along them
for iS = 1:size(sessions,1)
    monkey = sessions{iS,1};
    date = sessions{iS,2};
    
    s_idx = find(strcmpi(filedb.Date,date) & strcmpi(filedb.Monkey,monkey),1);
    
    if isempty(s_idx)% if it's not in the database yet, add it
        filedb = filedb_add(filedb);
        s_idx = find(strcmpi(filedb.Date,date) & strcmpi(filedb.Monkey,monkey),1);
    end
    
    epochs=filedb.Epochs{s_idx};
    filenums = cellfun(@(x) x(end-6:end-4),filedb.FileNames{s_idx},'UniformOutput',false);
    task=filedb.Task{s_idx};
    pert = filedb.Perturbation{s_idx};
    
    if str2double(date(1:4)) < 2016 % all Lab 3 data is 2015 or earlier
        lab = 3;
    else
        lab = 6;
    end
    
    % get array maps
    arrays = array_list{strcmpi(array_list(:,1),monkey) & cellfun(@(x) datenum(date) < datenum(x), array_list(:,2)),3};
    arrayMaps = array_list{strcmpi(array_list(:,1),monkey) & cellfun(@(x) datenum(date) < datenum(x), array_list(:,2)),4};
    
    if makeTD,
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
            cds.file2cds(fname,lab,['array' arrays{1}],['monkey' monkey],['task' task],'ranByMatt',['mapFile' arrayMaps{1}],'ignoreJumps')
            
            if length(arrays) > 1 % add in the PMd data
                fname = fullfile(rootDir,monkey,cerebusDataDir,date,[monkey '_' arrays{2} '_' task '_' pert '_' epoch '_' datestr(date,'mmddyyyy') '_' filenum '.nev']);
                cds.aliasList = {};
                % load spiking for PMd array
                cds.file2cds(fname,lab,['array' arrays{2}],['monkey' monkey],['task' task],'ranByMatt',['mapFile' arrayMaps{2}],'ignoreJumps')
            end
            
            % save it
            if ~exist(fullfile(rootDir,monkey,CDSDir,date),'dir')
                mkdir(fullfile(rootDir,monkey,CDSDir,date));
            end
            save(outname,'-v7.3','cds')
        else
            disp([epoch ' CDS already exists... loading...']);
            if makeTD
                load(outname,'cds');
            end
        end
        
        if makeTD
            meta.epoch = epoch;
            tdInputArgs.meta = meta;
            td = parseFileByTrial_cds(cds,tdInputArgs);
            
            trial_data = [trial_data, td];
        end
        clear td cds;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check quality of trials and prune some bad ones
    if makeTD
        bad_idx = cellfun(@(x) any(sqrt(x(:,1).^2+x(:,2).^2) > max_vel),{trial_data.vel});
        trial_data(bad_idx) = [];
        disp(['Pruning ' num2str(sum(bad_idx)) ' trials with crazy velocities...']);
        
        disp('Checking sorted units...');
        trial_data = getCommonUnits(trial_data);
        save(fullfile(rootDir,TDDir,[monkey '_' task '_' pert '_' date '.mat']),'trial_data');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add CDS info to filedb
    if somethingCDS % assume if this was false that there is nothing new to add
        filedb = filedb_add(filedb,fullfile(rootDir,monkey,CDSDir,date));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
