function trial_data = parseFileByTrial_bdf(data,inputArgs)
% inputArgs:
%   data: REQUIRED... my data struct
%   binSize: default 0.01 sec
%   dt: time resolution of kinematic signal in data, default 0.001 sec
%   extraTime: [time before, time after] beginning and end of trial (default [0.5 0.3] sec)
possibleArrays = {'M1','PMd'};

if isfield(inputArgs,'binSize'), binSize = inputArgs.binSize; else binSize = 0.01; end
if isfield(inputArgs,'dt'), dt = inputArgs.dt; else dt = 0.001; end
if isfield(inputArgs,'extraTime'), extraTime = inputArgs.extraTime; else extraTime = [0.5 0.3]; end

% see what array data is present
useArrays = intersect(possibleArrays,fieldnames(data));

% check if data file has the expected format of meta data
%   problem is that when I stitch files together some info is added
%   file_ids is added when multiple files are stitched together to identify which trials
if ~isfield(data.meta,'file_ids')
    data.meta.file_ids = ones(size(data.movement_table,1),1);
    data.meta.task = {data.meta.task};
    data.meta.perturbation = {data.meta.perturbation};
    data.meta.epoch = {data.meta.epoch};
end

mt = data.movement_table;

% loop along trials
trial_data = repmat(struct(),1,size(mt,1));
for iTrial = 1:size(mt,1)
    % add some meta data about the trial
    fileID = data.meta.file_ids(iTrial);
    trial_data(iTrial).monkey = data.meta.monkey;
    trial_data(iTrial).date = data.meta.recording_date;
    trial_data(iTrial).epoch = data.meta.epoch{fileID};
    trial_data(iTrial).task = data.meta.task{fileID};
    trial_data(iTrial).perturbation = data.meta.perturbation{fileID};
    trial_data(iTrial).target_direction = mt(iTrial,1);
    trial_data(iTrial).trialId = iTrial;
    
    switch lower(data.params.exp.angle_dir)
        case 'ccw'
            perturbation_direction = 1;
        case 'cw'
            perturbation_direction = -1;
    end
    
    % for VR: angle (negative is clockwise); for FF: [magnitude, direction] (negative is clockwise)
    switch lower(data.meta.perturbation{fileID})
        case 'ff'
            trial_data(iTrial).perturbation_info = [data.params.exp.force_magnitude, perturbation_direction*data.params.exp.force_angle];
        case 'vr'
            trial_data(iTrial).perturbation_info = perturbation_direction*data.params.exp.rotation_angle;
    end
    
    % find trial start/end times
    t_start = mt(iTrial,2) - extraTime(1);
    t_end = mt(iTrial,6) + extraTime(2);
    
    % get time vector for binned spikes
    t_bins = t_start:binSize:t_end;
    
    % put trial markers (target on etc) in bins for each spikes
    trial_data(iTrial).idx_target_on = find(histcounts(mt(iTrial,2),t_bins));
    trial_data(iTrial).idx_go_cue = find(histcounts(mt(iTrial,3),t_bins));
    trial_data(iTrial).idx_movement_on = find(histcounts(mt(iTrial,4),t_bins));
    trial_data(iTrial).idx_peak_speed = find(histcounts(mt(iTrial,5),t_bins));
    trial_data(iTrial).idx_reward = find(histcounts(mt(iTrial,6),t_bins));
    
    % get kinematics for that trial
    idx = data.cont.t >= t_start + binSize/2 & data.cont.t < t_end - binSize/2;
    trial_data(iTrial).pos = [decimate(data.cont.pos(idx,1),round(binSize/dt)) decimate(data.cont.pos(idx,2),round(binSize/dt))];
    trial_data(iTrial).vel = [decimate(data.cont.vel(idx,1),round(binSize/dt)) decimate(data.cont.vel(idx,2),round(binSize/dt))];
    trial_data(iTrial).acc = [decimate(data.cont.acc(idx,1),round(binSize/dt)) decimate(data.cont.acc(idx,2),round(binSize/dt))];
    if ~isempty(data.cont.force)
        trial_data(iTrial).force = [decimate(data.cont.force(idx,1),round(binSize/dt)) decimate(data.cont.force(idx,2),round(binSize/dt))];
    else
        trial_data(iTrial).force = [];
    end
    
    for iArray = 1:length(useArrays)
        binned_spikes = zeros(size(data.(useArrays{iArray}).sg,1),length(t_bins)-1);
        for unit = 1:size(data.(useArrays{iArray}).sg,1);
            % get the spikes for that cell in the current time window
            ts = data.(useArrays{iArray}).units(unit).ts;
            ts = ts(ts >= t_start & ts < t_end);
            
            binned_spikes(unit,:) = histcounts(ts,t_bins);
        end
        
        % save it as a sparse matrix
        trial_data(iTrial).([useArrays{iArray} '_spikes']) = sparse(binned_spikes);
        trial_data(iTrial).([useArrays{iArray} '_unit_guide']) = data.(useArrays{iArray}).sg;
    end
    clear binned_spikes;
end

