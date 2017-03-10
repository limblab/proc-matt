function filedb = filedb_add(filedb,cds_path)
% Create and maintain database table of sessions for Matt's analysis
%
% INPUTS
% (Pass in no inputs to initialize the filedb from dataSummary.m)
% filedb: (table) the file database
%   NOTE: Passing a filedb as only input will check that filedb against
%   dataSummary.m and initialize any missing sessions, but keep all other
%   info already present
% cds: (str) path to folder with the session CDS files to add CDS info for
%   a given session. Currently only supports one session at a time.
%
% OUTPUTS
% filedb: the filedb table with added info
%
% uses dataSummary.m to get master list of sessions
%
% Note: assumes all CDS files from a session will be together in a
% subfolder that does not include files from other sessions
%
% By default, the fileDB has the following entries:
% 'Monkey'
% 'Date'
% 'Perturbation'
% 'Task'
% 'Arrays'
% 'Direction'
% 'Notes'
% 'FFMagnitude'
% 'FFAngle'
% 'VRAngle'
% 'GRNumRewards'
%
% Once a CDS has been created, the following will be added
% 'Arrays': which arrays
% 'Units': count of neurons on that session
% 'FileNames': cell with filename of each file
% 'Reward': array with reward count for each file
% 'Abort': array with abort count for each file
% 'Failure': array with failure count for each file
% 'Incomplete': array with incomplete count for each file
% 'Epochs': cell array with names of epoch for each file
% 'Duration': array with length of each file
defineDirs;
%%%%%%%%%%%
% Some hard coded parameters
params.dbDir = dbDir;

if ~exist(fullfile(dbDir,'filedb.mat'),'file')
    filedb = table();
    save(fullfile(dbDir,'filedb.mat'),'filedb');
end

dataSummary;
params.cerebusDataDir = cerebusDataDir;
params.array_list = array_list;
params.exclude_units = [0 255];
if nargin < 1 % load session list from data summary and process
    disp('Initializing filedb...');
    did_something = true;
    
    filedb = table();
    for iFile = 1:size(sessionList,1)
        filedb = append_entry(filedb,sessionList(iFile,:),params);
    end
    
    % add root directory to use later
    filedb.Properties.UserData = dbDir;
    
elseif nargin == 1 % check against dataSummary.m for missing sessions and append them
    % check all monkeys and days
    idx = find(~(ismember(sessionList(:,1),filedb.Monkey) & ismember(sessionList(:,2),filedb.Date)));
    if ~isempty(idx) % session isn't added yet
        did_something = true;
        disp(['Appending ' num2str(length(idx)) ' missing sessions...']);
        for s_idx = idx'
            filedb = append_entry(filedb,sessionList(s_idx,:),params);
        end
    else
        did_something = false;
        disp('No missing sessions.');
    end
    
elseif nargin == 2 % add data from CDS
    did_something = true;
    
    d = dir(cds_path);
    filenames = {d(~[d.isdir]).name};
    
    % reorder filenames to fit chronological order
    filenames = filenames(cellfun(@(x) str2double(x(end-6:end-4)),filenames));
    
    [durations,rewards,aborts,failures,incompletes] = deal(zeros(1,length(filenames)));
    epochs = cell(1,length(filenames));
    for iFile = 1:length(filenames)
        % load CDS
        load(fullfile(cds_path,filenames{iFile}),'cds');
        
        % check if this session is already in the filedb
        s_idx = find(strcmpi(filedb.Date,datestr(cds.meta.dateTime,'yyyy-mm-dd')) & strcmpi(filedb.Monkey,cds.meta.monkey),1);
        if iFile == 1
            if isempty(s_idx) % if not, add from dataSummary
                disp('Session not present in filedb. Initializing...');
                dataSummary;
                s_idx = find(strcmpi(sessionList(:,1),cds.meta.monkey) & strcmpi(sessionList(:,2),datestr(cds.meta.dateTime,'yyyy-mm-dd')));
                
                filedb = append_entry(filedb,sessionList(s_idx,:),params);
            end
            
            arrays = strsplit(cds.meta.array,', ');
            
            units = zeros(1,length(arrays));
            for iArray = 1:length(arrays)
                units(iArray) = sum(~ismember([cds.units.ID],params.exclude_units) & strcmpi({cds.units.array},arrays{iArray}));
            end
        else
            % check number of units against original to ensure
            for iArray = 1:length(arrays)
                temp_units = sum(~ismember([cds.units.ID],params.exclude_units) & strcmpi({cds.units.array},arrays{iArray}));
                if temp_units ~= units(iArray), warning('Unit counts are different between files...'); end
                units(iArray) = min([units(iArray),temp_units]);
            end
        end
        
        epochs{iFile} = filenames{iFile}(end-18:end-17);
        
        durations(iFile) = cds.meta.duration;
        rewards(iFile) = cds.meta.numReward;
        aborts(iFile) = cds.meta.numAbort;
        failures(iFile) = cds.meta.numFail;
        incompletes(iFile) = cds.meta.numIncomplete;
    end
    
    % now add to filedb
    filedb.FileNames{s_idx} = filenames;
    filedb.Epochs{s_idx} = epochs;
    filedb.Duration{s_idx} = durations;
    filedb.Units{s_idx} = units;
    filedb.Reward{s_idx} = rewards;
    filedb.Abort{s_idx} = aborts;
    filedb.Failure{s_idx} = failures;
    filedb.Incomplete{s_idx} = incompletes;
else
    error('Too many inputs.');
end

if did_something
    disp('Saving...')
    save(fullfile(dbDir,'filedb.mat'),'filedb');
    disp('Done.');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add a row to the database table
function filedb = append_entry(filedb,s,params)

dbDir = params.dbDir;
cerebusDataDir = params.cerebusDataDir;
array_list = params.array_list;

switch lower(s{3})
    case 'ff'
        ff_mag = s{7};
        ff_ang = s{6};
        vr_ang = NaN;
        gr_rewards = NaN;
    case 'vr'
        ff_mag = NaN;
        ff_ang = NaN;
        vr_ang = s{6};
        gr_rewards = NaN;
    case 'gr'
        ff_mag = NaN;
        ff_ang = NaN;
        vr_ang = s{6};
        gr_rewards = s{7};
    case 'cs'
        ff_mag = NaN;
        ff_ang = NaN;
        vr_ang = NaN;
        gr_rewards = NaN;
end

% get filenames for all associated files
d = dir(fullfile(dbDir,s{1},cerebusDataDir,s{2}));
d=d(~[d.isdir]);
fn = {d(~[d.isdir]).name};
% get end of filenames (we want '001', '002', etc)
fn = fn(strcmpi(cellfun(@(x) x(end-2:end),fn,'UniformOutput',false),'ccf'));
e = unique(cellfun(@(x) x(end-6:end-4),fn,'UniformOutput',false));

% loop along the unique numbers and get the associated epoch
[filenames,epochs] = deal(cell(1,length(e)));
for i = 1:length(e)
    % there will be multiple files with the same suffix if there
    % are dual array recordings so just take the first one
    idx = find(strcmpi(cellfun(@(x) x(end-6:end-4),fn,'UniformOutput',false),e{i}));
    epochs{i} = fn{idx(1)}(end-18:end-17);
    filenames{i} = [s{1} '_' s{4} '_' s{3} '_' datestr(s{2},'mmddyyyy') '_' e{i} '.mat'];
end

% now get arrays for this monkey from array_list
arrays = array_list{strcmpi(array_list(:,1),s{1}) & cellfun(@(x) datenum(s{2}) < datenum(x), array_list(:,2)),3};

session_info = table( ...
    s(1), s(2), s(3), s(4), ...
    {arrays}, {filenames}, {epochs}, cell(1), ...
    s(5), ff_mag, ff_ang, vr_ang, gr_rewards, ...
    cell(1),cell(1),cell(1),cell(1),cell(1),s(8), ...
    'VariableNames', {'Monkey','Date','Perturbation','Task', ...
    'Arrays','FileNames','Epochs','Duration', ...
    'Direction','FFMagnitude','FFAngle','VRAngle','GRNumRewards', ...
    'Units','Reward','Abort','Failure','Incomplete','Notes'});

filedb = vertcat(filedb,session_info);
end