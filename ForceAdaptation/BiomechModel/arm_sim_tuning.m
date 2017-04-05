%% Do onset to peak tuning
coordinate_frame = 'movement';
% blocks = {'BL',[0 1]; 'AD',[0 0.33]; 'AD',[0.33 0.66]; 'AD',[0.66 1]; 'WO',[0 0.33]; 'WO',[0.33 0.66]; 'WO',[0.66 1]};
blocks = {'BL',[0 1]; 'AD',[0.66 1]; 'WO',[0.66 1]};
% blocks = {'BL',[0 1]};
angle_bin_size = 45*pi/180;

if tune_real_neurons
    disp('USING REAL NEURONS!');
    n = size(trial_data(1).M1_spikes,2);
else
    n = params.num_neurons;
end

tc_params.coordinate_frame = coordinate_frame;
tc_params.blocks = blocks;
tc_params.angle_bin_size = angle_bin_size;
params.tc_params = tc_params;

theta = zeros(length(sim_data),1);
for iTrial = 1:length(sim_data)
    switch lower(coordinate_frame)
        case 'movement'
            if ~isempty(trial_data(iTrial).idx_peak_speed) && (trial_data(iTrial).idx_peak_speed - trial_data(iTrial).idx_movement_on) >= 10
                idx = trial_data(iTrial).idx_movement_on:trial_data(iTrial).idx_peak_speed;
                theta(iTrial) = atan2(trial_data(iTrial).pos(idx(end),2)-trial_data(iTrial).pos(idx(1),2),trial_data(iTrial).pos(idx(end),1)-trial_data(iTrial).pos(idx(1),1));
            else
                theta(iTrial) = NaN;
            end
        case 'target'
            theta(iTrial) = trial_data(iTrial).target_direction;
    end
end
bad_trials = isnan(theta)';
theta = bin_angles(theta,angle_bin_size);

clear tc_data;

for i = 1:length(use_models)
    
    fr = zeros(length(sim_data),n);
    ma = zeros(length(sim_data),params.how_many_muscles);
    for iTrial = 1:length(sim_data)
        if ~isempty(trial_data(iTrial).idx_peak_speed) && (trial_data(iTrial).idx_peak_speed - trial_data(iTrial).idx_movement_on) >= 10
            idx = trial_data(iTrial).idx_movement_on:trial_data(iTrial).idx_peak_speed;
            
            idx = idx - mean_lag;
            
            if tune_real_neurons
                fr(iTrial,:) = sum(trial_data(iTrial).M1_spikes(idx,:),2)./(dt*length(idx));
            else
                fr(iTrial,:) = sum(sim_data(iTrial).([use_models{i} '_neurons'])(idx,:),1)./(dt*length(idx));
            end
            
            ma(iTrial,:) = mean(sim_data(iTrial).muscles(idx,:),1);
        end
    end
    
    for j = 1:size(blocks,1)
        idx = find(strcmpi({trial_data.epoch},blocks{j,1}) & ~bad_trials);
        idx = idx(floor(blocks{j,2}(1)*length(idx))+1:floor(blocks{j,2}(2)*length(idx)));
        
        if j == 1 % do anova of direction in baseline
            p = zeros(1,size(fr,2));
            for k = 1:size(fr,2)
                p(k) = anovan(fr(idx,k),theta(idx),'display','off');
            end
        end
        
        [tunCurves,confBounds,rs,boot_pds,boot_mds,boot_bos] = regressTuningCurves(fr(idx,:),theta(idx),{'bootstrap',1000,0.95},'doplots',false,'doparallel',true,'domeanfr',true);
        tc_data(j).(use_models{i}).tc = tunCurves;
        tc_data(j).(use_models{i}).cb = confBounds;
        tc_data(j).(use_models{i}).rs = rs;
        tc_data(j).(use_models{i}).boot_pds = boot_pds;
        tc_data(j).(use_models{i}).boot_mds = boot_mds;
        tc_data(j).(use_models{i}).boot_bos = boot_bos;
        tc_data(j).(use_models{i}).anova_p = p;
    end
    
end

%% Do sliding window tuning
if do_sliding_window_bullshit
    coordinate_frame = 'movement';
    blocks = {'BL',[0 1]; 'AD',[0.33 1]};
    angle_bin_size = 45*pi/180;
    time_delay = 0;
    hold_time = 50;
    divide_time = [0.3 0.05];
    
    numBlocks = floor(( 1 + divide_time(2) - divide_time(1) ) / divide_time(2));
    
    sw_params.coordinate_frame = coordinate_frame;
    sw_params.blocks = blocks;
    sw_params.angle_bin_size = angle_bin_size;
    sw_params.time_delay = time_delay;
    sw_params.hold_time = hold_time;
    sw_params.divide_time = divide_time;
    
    params.sw_params = sw_params;
    
    %
    disp('HEY YOURE ONLY DOING SIX BLOCKS');
    clear sw_data;
    for iBlock = 1:6%numBlocks
        clear temp temp2;
        for i = 1:length(use_models)
            theta = zeros(length(sim_data),1);
            
            fr = zeros(length(sim_data),n);
            f = zeros(length(sim_data),1);
            for iTrial = 1:length(sim_data)
                if ~isempty(trial_data(iTrial).idx_peak_speed) && (trial_data(iTrial).idx_peak_speed - trial_data(iTrial).idx_movement_on) >= 10
                    
                    tstart = trial_data(iTrial).idx_movement_on-time_delay; % time that reach starts
                    tdur = ( trial_data(iTrial).idx_trial_end - hold_time + time_delay ) - tstart; % time duration of reach
                    idx = floor(tstart + (iBlock-1)*divide_time(2)*tdur) : floor(tstart + (divide_time(1)+(iBlock-1)*divide_time(2)) * tdur);
                    
                    try % SHORT TERM HACK FIX THIS
                        %theta(iTrial) = atan2(trial_data(iTrial).vel(idx(end),2),trial_data(iTrial).vel(idx(end),1));
                        theta(iTrial) = atan2(trial_data(iTrial).pos(idx(end),2)-trial_data(iTrial).pos(idx(1),2),trial_data(iTrial).pos(idx(end),1)-trial_data(iTrial).pos(idx(1),1));
                        f(iTrial) = rms(hypot(trial_data(iTrial).force(idx,1),trial_data(iTrial).force(idx,2)));
                        
                        idx = idx - mean_lag;
                        
                        if tune_real_neurons
                            fr(iTrial,:) = sum(trial_data(iTrial).M1_spikes(idx,:),2)./(dt*length(idx));
                        else
                            fr(iTrial,:) = sum(sim_data(iTrial).([use_models{i} '_neurons'])(idx,:),1)./(dt*length(idx));
                        end
                    catch
                        disp('FUCK ERRORS');
                        theta(iTrial) = NaN;
                    end
                else
                    theta(iTrial) = NaN;
                end
            end
            
            bad_trials = isnan(theta)';
            
            theta = bin_angles(theta,angle_bin_size);
            
            for j = 1:length(blocks)
                idx = find(strcmpi({trial_data.epoch},blocks{j,1}) & ~bad_trials);
                idx = idx(floor(blocks{j,2}(1)*length(idx))+1:floor(blocks{j,2}(2)*length(idx)));
                
                [tunCurves,confBounds,rs,boot_pds,boot_mds,boot_bos] = regressTuningCurves(fr(idx,:),theta(idx),{'bootstrap',1000,0.95},'doplots',false,'doparallel',true,'domeanfr',true);
                temp(j).(use_models{i}).tc = tunCurves;
                temp(j).(use_models{i}).cb = confBounds;
                temp(j).(use_models{i}).rs = rs;
                temp(j).(use_models{i}).boot_pds = boot_pds;
                temp(j).(use_models{i}).boot_mds = boot_mds;
                temp(j).(use_models{i}).boot_bos = boot_bos;
                
                temp2(j).f = f(idx);
                temp2(j).theta = theta(idx);
                temp2(j).bad_trials = bad_trials(idx);
            end
        end
        sw_data(iBlock).tc_data = temp;
        sw_data(iBlock).data = temp2;
    end
else
    disp('SKIPPING SLIDING WINDOW THING!!!!!!');
    sw_data = [];
end


