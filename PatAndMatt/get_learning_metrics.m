function metric = get_learning_metrics(trial_data, which_metric, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% will compute and plot behavioral adaptation metrics
%
% INPUTS:
%   trial_data: the struct
%   metric: which metric
%           'angle': angular takeoff error
%           'corr': speed profile correlation coefficients
%           'time': time to target
%   params: struct with the following options
%           'result_codes': which to include. Default is 'R', 'I'.
%           'corr_samples': how many datapoints to interpolate trajectory onto for corr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_codes = {'R','I'};
corr_samples = 1000;
if nargin == 3
    if isfield(params,'result_codes'), result_codes = params.result_codes; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
velorpos = 'pos';

if isfield(trial_data(1),'result')
    trial_data = trial_data(ismember({trial_data.result},result_codes));
else
    disp('Result field not found...');
end

utheta = unique([trial_data.target_direction]);
metric = zeros(length(trial_data),1);
switch lower(which_metric)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'angle'
        
        % get baseline error to each target
        bl_metric = zeros(length(utheta),1);
        for iDir = 1:length(utheta)
            bl_idx = find(getTDidx(trial_data,'epoch','bl','target_direction',utheta(iDir)));
            temp = zeros(length(bl_idx),1);
            for iTrial = 1:length(bl_idx)
                switch velorpos
                    case 'vel'
                        temp(iTrial) = angleDiff(trial_data(bl_idx(iTrial)).target_direction, ...
                            atan2(trial_data(bl_idx(iTrial)).vel(trial_data(bl_idx(iTrial)).idx_peak_speed,2), ...
                            trial_data(bl_idx(iTrial)).vel(trial_data(bl_idx(iTrial)).idx_peak_speed,1)), ...
                            true,true);
                    case 'pos'
                        temp(iTrial) = angleDiff(trial_data(bl_idx(iTrial)).target_direction, ...
                            atan2(trial_data(bl_idx(iTrial)).vel(trial_data(bl_idx(iTrial)).idx_peak_speed,2) - ...
                            trial_data(bl_idx(iTrial)).vel(trial_data(bl_idx(iTrial)).idx_movement_on,2), ...
                            trial_data(bl_idx(iTrial)).vel(trial_data(bl_idx(iTrial)).idx_peak_speed,1) - ...
                            trial_data(bl_idx(iTrial)).vel(trial_data(bl_idx(iTrial)).idx_movement_on,1)), ...
                            true,true);
                end
            end
            bl_metric(iDir) = circular_mean(temp);
        end, clear temp;
        
        % get velocity at time of peak speed
        for iTrial = 1:length(trial_data)
            switch velorpos
                case 'vel'
                    temp = angleDiff(trial_data(iTrial).target_direction, ...
                        atan2(trial_data(iTrial).vel(trial_data(iTrial).idx_peak_speed,2), ...
                        trial_data(iTrial).vel(trial_data(iTrial).idx_peak_speed,1)), ...
                        true,true);
                case 'pos'
                    temp = angleDiff(trial_data(iTrial).target_direction, ...
                        atan2(trial_data(iTrial).vel(trial_data(iTrial).idx_peak_speed,2) - ...
                        trial_data(iTrial).vel(trial_data(iTrial).idx_movement_on,2), ...
                        trial_data(iTrial).vel(trial_data(iTrial).idx_peak_speed,1) - ...
                        trial_data(iTrial).vel(trial_data(iTrial).idx_movement_on,1)), ...
                        true,true);
            end
            iDir = utheta==trial_data(iTrial).target_direction;
            metric(iTrial) = angleDiff(bl_metric(iDir),temp,true,true);
        end, clear temp;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'corr'
        % get baseline trace to each target
        bl_metric = zeros(length(utheta),corr_samples,2);
        for iDir = 1:length(utheta)
            bl_idx = find(getTDidx(trial_data,'epoch','bl','target_direction',utheta(iDir)));
            bl_temp = zeros(length(bl_idx),2,corr_samples);
            for iTrial = 1:length(bl_idx)
                idx = trial_data(bl_idx(iTrial)).idx_movement_on-10:trial_data(bl_idx(iTrial)).idx_trial_end-30;
                temp = trial_data(bl_idx(iTrial)).vel;
                %                 temp(:,1) = smoothSpikesForPCA(temp(:,1),1,10);
                %                 temp(:,2) = smoothSpikesForPCA(temp(:,2),1,10);
                bl_temp(iTrial,1,:) = interp1(1:length(idx),temp(idx,1),linspace(1,length(idx),corr_samples));
                bl_temp(iTrial,2,:) = interp1(1:length(idx),temp(idx,2),linspace(1,length(idx),corr_samples));
            end
            bl_metric(iDir,:,:) = squeeze(mean(bl_temp,1))';
        end, clear temp bl_temp;
        
        for iTrial = 1:length(trial_data)
            idx = trial_data(iTrial).idx_movement_on-10:trial_data(iTrial).idx_trial_end-30;
            iDir = utheta==trial_data(iTrial).target_direction;
            
            temp = trial_data(iTrial).vel;
            if strcmpi(trial_data(iTrial).epoch,'ad')
                temp(:,1) = smoothSpikesForPCA(temp(:,1),1,10);
                temp(:,2) = smoothSpikesForPCA(temp(:,2),1,10);
            end
            temp = [interp1(1:length(idx),temp(idx,1),linspace(1,length(idx),corr_samples))', ...
                interp1(1:length(idx),temp(idx,2),linspace(1,length(idx),corr_samples))'];
            metric(iTrial) = corr2(squeeze(bl_metric(iDir,:,:)),temp)^2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'time'
        error('Time to target not implemented.');
        
    otherwise
        error('metric not recognized.');
end


