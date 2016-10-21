function trial_data = parseFileByTrial_cds(cds,inputArgs)
%   data: REQUIRED... cell array of CDS objects. Will concatenate together.
%
% inputArgs:
%   meta: .perturbation, .epoch, .angle_dir, .rotation_angle, .force_magnitude, .force_angle (REQUIRED)
%   trialResults: which reward codes to use ('R','A','F','I')
%   binSize: default 0.01 sec
%   extraTime: [time before, time after] beginning and end of trial (default [0.5 0.3] sec)
if isfield(inputArgs,'trialResults'), trialResults = inputArgs.trialResults; else trialResults = {'R'}; end
if isfield(inputArgs,'excludeUnits'), excludeUnits = inputArgs.excludeUnits; else excludeUnits = [0,255]; end
if isfield(inputArgs,'binSize'), binSize = inputArgs.binSize; else binSize = 0.01; end
if isfield(inputArgs,'extraTime'), extraTime = inputArgs.extraTime; else extraTime = [0.1 0.1]; end

% some hard coded parameters
min_ds = 0.3; % minimum diff(speed) to find movement onset
min_trial_time = 0.1; % minimum time from start to end

% do some input processing
if ~iscell(trialResults), trialResults = {trialResults}; end

% see what array data is present
arrays = strsplit(cds.meta.array,', ');

% get info on neurons
unit_idx = cell(1,length(arrays));
for iArray = 1:length(arrays)
    unit_idx{iArray} = find(~ismember([cds.units.ID],excludeUnits) & strcmpi({cds.units.array},arrays{iArray}));
end

idx_trials = find(ismember(cds.trials.result,trialResults));

dt_kin = cds.kin.t(2)-cds.kin.t(1);
dt_force = cds.force.t(2)-cds.force.t(1);

% loop along trials
bad_idx = [];
trial_data = repmat(struct(),1,length(idx_trials));
for i = 1:length(idx_trials)
    iTrial = idx_trials(i);
    
    if cds.trials.endTime(iTrial)-cds.trials.startTime(iTrial) > min_trial_time
        
        % add some meta data about the trial
        trial_data(i).monkey = cds.meta.monkey;
        trial_data(i).date = datestr(cds.meta.dateTime,'mm-dd-yyyy');
        trial_data(i).epoch = inputArgs.meta.epoch;
        trial_data(i).task = cds.meta.task;
        trial_data(i).perturbation = inputArgs.meta.perturbation;
        trial_data(i).target_direction = cds.trials.tgtDir(iTrial);
        trial_data(i).trial_id = iTrial;
        trial_data(i).result = cds.trials.result(iTrial);
        
        switch lower(inputArgs.meta.angle_dir)
            case 'ccw'
                perturbation_direction = 1;
            case 'cw'
                perturbation_direction = -1;
        end
        
        % for VR: angle (negative is clockwise); for FF: [magnitude, direction] (negative is clockwise)
        switch lower(inputArgs.meta.perturbation)
            case 'ff'
                trial_data(i).perturbation_info = [inputArgs.meta.force_magnitude, perturbation_direction*inputArgs.meta.force_angle];
            case 'vr'
                trial_data(i).perturbation_info = perturbation_direction*inputArgs.meta.rotation_angle;
        end
        
        % find trial start/end times
        t_start = cds.trials.startTime(iTrial) - extraTime(1);
        t_end = cds.trials.endTime(iTrial) + extraTime(2);
        
        % get time vector for binned spikes
        t_bins = t_start:binSize:t_end;
        
        % get kinematics for that trial
        idx = cds.kin.t >= t_start + binSize/2 & cds.kin.t < t_end - binSize/2;
        trial_data(i).pos = [decimate(cds.kin.x(idx),round(binSize/dt_kin)) decimate(cds.kin.y(idx),round(binSize/dt_kin))];
        trial_data(i).vel = [decimate(cds.kin.vx(idx),round(binSize/dt_kin)) decimate(cds.kin.vy(idx),round(binSize/dt_kin))];
        trial_data(i).acc = [decimate(cds.kin.ax(idx),round(binSize/dt_kin)) decimate(cds.kin.ay(idx),round(binSize/dt_kin))];
        if ~isempty(cds.force)
            trial_data(i).force = [decimate(cds.force.fx(idx),round(binSize/dt_force)) decimate(cds.force.fy(idx),round(binSize/dt_force))];
        else
            trial_data(i).force = [];
        end
        
        
        % put trial markers (target on etc) in bins for each spikes
        if ~isnan(cds.trials.startTime(iTrial)) && ~isnan(cds.trials.goCueTime(iTrial))
            trial_data(i).idx_trial_start = find(histcounts(cds.trials.startTime(iTrial),t_bins));
            trial_data(i).idx_target_on = find(histcounts(cds.trials.tgtOnTime(iTrial),t_bins));
            trial_data(i).idx_go_cue = find(histcounts(cds.trials.goCueTime(iTrial),t_bins));
            
            % use velocity to find bin corresponding to movement onset, movement offset, and peak speed
            s = sqrt(trial_data(i).vel(:,1).^2 + trial_data(i).vel(:,2).^2);
            ds = [0; diff(s)];
            dds = [0; diff(ds)];
            peaks = [dds(1:end-1)>0 & dds(2:end)<0; 0];
            mvt_peak = find(peaks & (1:length(peaks))' > trial_data(i).idx_go_cue & ds > min_ds, 1, 'first');
            
            if ~isempty(mvt_peak)
                thresh = ds(mvt_peak)/2;                             % Threshold is half max of acceleration peak
                on_idx = find(ds<thresh & (1:length(ds))'<mvt_peak,1,'last');
                % find movement peak as maximum velocity
                s(1:on_idx) = 0;
                [~, peak_idx] = max(s);
                
                trial_data(i).idx_movement_on = on_idx;
                trial_data(i).idx_peak_speed = peak_idx;
                %     trial_data(i).idx_movement_off = off_idx;
                trial_data(i).idx_trial_end = find(histcounts(cds.trials.endTime(iTrial),t_bins));
                
                if trial_data(i).idx_peak_speed < trial_data(i).idx_trial_end
                    for iArray = 1:length(arrays)
                        binned_spikes = zeros(size(unit_idx{iArray},1),length(t_bins)-1);
                        sg = zeros(length(unit_idx{iArray}),2);
                        for unit = 1:length(unit_idx{iArray})
                            % get the spikes for that cell in the current time window
                            ts = cds.units(unit_idx{iArray}(unit)).spikes.ts;
                            ts = ts(ts >= t_start & ts < t_end);
                            
                            binned_spikes(unit,:) = histcounts(ts,t_bins);
                            sg(unit,:) = [cds.units(unit_idx{iArray}(unit)).chan, cds.units(unit_idx{iArray}(unit)).ID];
                        end
                        
                        %   transpose binned_spikes to be consistent with kin
                        trial_data(i).([arrays{iArray} '_spikes']) = binned_spikes';
                        trial_data(i).([arrays{iArray} '_unit_guide']) = sg;
                    end
                    clear binned_spikes;
                else
                    bad_idx = [bad_idx,i];
                end
            else
                bad_idx = [bad_idx,i];
            end
        else
            bad_idx = [bad_idx, i];
        end
    else
        bad_idx = [bad_idx, i];
    end
end

disp(['Pruning ' num2str(length(bad_idx)) ' trials with bad trial info...']);
trial_data(bad_idx) = [];
