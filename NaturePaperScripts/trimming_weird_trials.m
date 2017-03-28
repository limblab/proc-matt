disp('HEY! Im doing that weird trial trimming thing to make my new data struct like my old data structs.');

min_ds = 0.3; % minimum diff(speed) to find movement onset

bad_idx = find(cellfun(@(x) any(sqrt(x(:,1).^2+x(:,2).^2) > 50),{trial_data.vel}));

for i = 1:length(trial_data)
    if trial_data(i).idx_trial_end - trial_data(i).idx_trial_start > 10
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
            
            if isempty(peak_idx) || isempty(on_idx) || isnan(peak_idx) || isnan(on_idx)
                bad_idx = [bad_idx,i];
            else
                if trial_data(i).idx_peak_speed >= trial_data(i).idx_trial_end
                    bad_idx = [bad_idx, i];
                end
            end
        else
            bad_idx = [bad_idx, i];
        end
    else
        bad_idx = [bad_idx,i];
    end
end

bad_idx = unique(bad_idx);

trial_data(bad_idx) = [];