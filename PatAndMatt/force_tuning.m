% This will fit bootstrapped PDs to both hand direction (kinematic
% velocity) and computed "hand force" (a la Vaadia)
close all;
clc;

use_array = 'M1';
% how to partition trials
%  First element of each row is the epoch
%  Second element is a vector specifying proportional [start,end] of trials
tuning_blocks = {'BL',[0 1]; ...
    'AD',[0 0.33]; ...
    'AD',[0.33 0.66]; ...
    'AD',[0.66 1]; ...
    'WO',[0 0.5]; ...
    'WO',[0.5 1]};

% some parameters for bootstrapping
num_iterations = 1000;
conf_level = 0.95;
dt = 0.01; % binned data bin size in sec

% total force = mass * acc - force
m = 0.7; %kg
idx_shift = 10; % how many 10ms bins to shift neural data

doPos = true;
doVel = true;
doMP = true;
nr_operation = 'replace'; % 'add' or 'replace'

%% Get information about the available data
dates = unique({trial_data.date});
exclude_dates = {'2015-06-29'};
dates = setdiff(dates,exclude_dates);

%% Get directions and firing rates in the window
all_d = [];
all_sig = [];
all_rs = [];

if ~exist('nr','var')
    nr = [];
    disp('Initializing neural results struct...');
end
for iFile = 1:length(dates)
    clc;
    disp(['File ' num2str(iFile) ' of ' num2str(length(dates))]);
    date_trial_inds = find(get_trial_data_indices(trial_data,'date',dates{iFile}));
    
    fr = zeros(length(date_trial_inds),size(trial_data(date_trial_inds(1)).([use_array '_spikes']),1));
    theta_p = zeros(length(date_trial_inds),1);
    theta_v = zeros(length(date_trial_inds),1);
    theta_f = zeros(length(date_trial_inds),1);
    target_directions = zeros(length(date_trial_inds),1);
    for iInd = 1:length(date_trial_inds)
        
        iTrial = date_trial_inds(iInd);
        
        target_directions(iTrial) = trial_data(iTrial).target_direction;
        if isempty(trial_data(iTrial).idx_peak_speed) || (trial_data(iTrial).idx_peak_speed <= trial_data(iTrial).idx_movement_on)
            [theta_p(iInd),theta_v(iInd),theta_f(iInd)] = deal(NaN);
            fr(iInd,:) = NaN(1,size(fr,2));
        else
            % determine force constant for this trial
            if strcmpi(trial_data(iTrial).epoch,'AD')
                k = trial_data(iTrial).perturbation_info(1); %Ns/cm
                th = trial_data(iTrial).perturbation_info(2); %rad
                if isnan(th)
                    error('NOOOOOO!');
                end
            else
                k = 0;
                th = 0;
            end
            
            % get index between movement onset and peak speed
            %         idx = trial_data(iTrial).idx_movement_on:trial_data(iTrial).idx_movement_on+15;
            idx = trial_data(iTrial).idx_movement_on:trial_data(iTrial).idx_peak_speed;
            
            % get velocity vector
            p = trial_data(iTrial).pos(idx,:);
            v = trial_data(iTrial).vel(idx,:);
            %a = trial_data(iTrial).acc(idx,:);
            a = [gradient(v(:,1),1/100),gradient(v(:,2),1/100)];
            f = trial_data(iTrial).force(idx,:);
            
            % get angle of hand movement
            theta_p(iInd) = atan2(p(end,2)-p(1,2),p(end,1)-p(1,1));
            theta_v(iInd) = atan2(v(end,2),v(end,1));
            
            if strcmpi(trial_data(iTrial).epoch,'AD')
                % get robot force vector
                f_robot = zeros(size(v,1),2);
                f_hand = zeros(size(v,1),2);
                f_plan = zeros(size(v,1),2);
                for j = 1:size(v,1)
                    f_robot(j,:) = k.*[cos(th),-sin(th);sin(th),cos(th)]*v(j,:)';
                    %f_hand(j,:) = (1/100).*m.*a(j,:); % 1/100 converts cm/s^2 to m/s^2
                    f_hand(j,:) = f(j,:);
                    f_plan(j,:) = f_hand(j,:)' - f_robot(j,:)';
                end
                
                % find time bin corresponding to peak of acceleration
                [~,idx_peak] = max( hypot(f_plan(:,1),f_plan(:,2)) );
                theta_f(iInd) = atan2(f_plan(idx_peak,2),f_plan(idx_peak,1));
            else
                theta_f(iInd) = theta_v(iInd);
            end
            
            % get firing rates in same window for each cell
            spike_counts = full(trial_data(iTrial).([use_array '_spikes']));
            fr(iInd,:) = sum(spike_counts(:,idx-idx_shift),2)/(length(idx)*dt);
        end
    end
    
    % get list of epochs for filtering
    epochs = {trial_data.epoch};
    % filter out bad trials
    bad_trials = isnan(theta_p);
    theta_p = theta_p(~bad_trials);
    theta_v = theta_v(~bad_trials);
    theta_f = theta_f(~bad_trials);
    fr = fr(~bad_trials,:);
    target_directions = target_directions(~bad_trials);
    epochs = epochs(~bad_trials);
    
    % Fit cosine tuning curves to predicted force vector and actual hand velocity
    [rs_p,pds_p,mds_p,bos_p,rs_v,rs_f,pds_v,mds_v,bos_v,pds_f,mds_f,bos_f] = deal(cell(1,size(tuning_blocks,1)));
    for iBlock = 1:size(tuning_blocks,1)
        disp(['Block ' num2str(iBlock) ' of ' num2str(size(tuning_blocks,1)) '...']);
        idx = find(strcmpi(epochs,tuning_blocks{iBlock,1}));
        idx = idx( ceil(tuning_blocks{iBlock,2}(1)*length(idx))+1:ceil(tuning_blocks{iBlock,2}(2)*length(idx)) );
        
        if doPos
            [~,~,rs_p,pds_p,mds_p,bos_p] = regressTuningCurves(fr(idx,:),theta_p(idx),{'bootstrap',num_iterations,conf_level},'doParallel',true);
        end
        if doVel
            [~,~,rs_v,pds_v,mds_v,bos_v] = regressTuningCurves(fr(idx,:),theta_v(idx),{'bootstrap',num_iterations,conf_level},'doParallel',true);
        end
        if doMP
            [~,~,rs_f,pds_f,mds_f,bos_f] = regressTuningCurves(fr(idx,:),theta_f(idx),{'bootstrap',num_iterations,conf_level},'doParallel',true);
        end
        
        for unit = 1:size(fr,2)
            % Add this cosine fit to my neuron database
            analysis_info.date = dates{iFile};
            analysis_info.monkey = trial_data(date_trial_inds(1)).monkey;
            analysis_info.task = trial_data(date_trial_inds(1)).task;
            analysis_info.perturbation = trial_data(date_trial_inds(1)).perturbation;
            analysis_info.array = use_array;
            analysis_info.neuron_id = unit;
            analysis_info.analysis_type = 'tuning';
            analysis_info.method = 'regression';
            
            if doMP
                analysis_info.name = ['onpeak_mp_' tuning_blocks{iBlock,1} '_' num2str(tuning_blocks{iBlock,2}(1)) '_' num2str(tuning_blocks{iBlock,2}(2))];
                analysis_info.notes = ['Block ' tuning_blocks{iBlock,1} ', Trials ' num2str(tuning_blocks{iBlock,2}(1)) ':' num2str(tuning_blocks{iBlock,2}(2)) '; motor plan coordinates; movement onset to peak speed window'];
                analysis_info.bootstrap = [bos_f(unit,:); mds_f(unit,:); pds_f(unit,:)];
                analysis_info.r_squared = rs_f(unit,:);
                analysis_info.fr = fr(idx,unit);
                analysis_info.theta = theta_f(idx);
                nr = neuron_results_struct(nr_operation,nr,analysis_info);
            end
            
            if doPos
                analysis_info.name = ['onpeak_pos_' tuning_blocks{iBlock,1} '_' num2str(tuning_blocks{iBlock,2}(1)) '_' num2str(tuning_blocks{iBlock,2}(2))];
                analysis_info.notes = ['Block ' tuning_blocks{iBlock,1} ', Trials ' num2str(tuning_blocks{iBlock,2}(1)) ':' num2str(tuning_blocks{iBlock,2}(2)) '; position coordinates; movement onset to peak speed window'];
                analysis_info.bootstrap = [bos_p(unit,:); mds_p(unit,:); pds_p(unit,:)];
                analysis_info.r_squared = rs_p(unit,:);
                analysis_info.theta = theta_p(idx);
                nr = neuron_results_struct(nr_operation,nr,analysis_info);
            end
            
            if doVel
                analysis_info.name = ['onpeak_vel_' tuning_blocks{iBlock,1} '_' num2str(tuning_blocks{iBlock,2}(1)) '_' num2str(tuning_blocks{iBlock,2}(2))];
                analysis_info.notes = ['Block ' tuning_blocks{iBlock,1} ', Trials ' num2str(tuning_blocks{iBlock,2}(1)) ':' num2str(tuning_blocks{iBlock,2}(2)) '; velocity coordinates; movement onset to peak speed window'];
                analysis_info.bootstrap = [bos_v(unit,:); mds_v(unit,:); pds_v(unit,:)];
                analysis_info.r_squared = rs_v(unit,:);
                analysis_info.theta = theta_v(idx);
                nr = neuron_results_struct(nr_operation,nr,analysis_info);
            end
        end
    end
end

save('F:\trial_data_files\CO_FF_neuron_results.mat','nr')
clc;
disp('Done.');

%%
cell_idx = find(strcmpi({nr.monkey},'Mihili') | strcmpi({nr.monkey},'Chewie'));

tw = 'onpeak';
use_names = {'pos','mp'};
plot_colors = {'k','r','b'};

r2_min = 0.4;

x_bins = -180:10:180;

% first check if cells are tuned
is_tuned = zeros(size(tuning_blocks,1),length(cell_idx),length(use_names));
for i = 1:length(cell_idx)
        unit = cell_idx(i);
    
    for iName = 1:length(use_names)
        n = use_names{iName};
        
        idx_tune_bl = strcmpi({nr(unit).tuning.name},[tw '_' n '_' tuning_blocks{1,1} '_' num2str(tuning_blocks{1,2}(1)) '_' num2str(tuning_blocks{1,2}(2))]);
        
        for iBlock = 1:size(tuning_blocks,1)
            
            % get velocity indices
            idx_tune = strcmpi({nr(unit).tuning.name},[tw '_' n '_' tuning_blocks{iBlock,1} '_' num2str(tuning_blocks{iBlock,2}(1)) '_' num2str(tuning_blocks{iBlock,2}(2))]);
            
            is_tuned(iBlock,i,iName) = mean(nr(unit).tuning(idx_tune).r_squared) > r2_min;
        end
    end
end

use_tuning_blocks = 1:6;
always_tuned = zeros(length(cell_idx),length(use_tuning_blocks));
for iBlock = 1:length(use_tuning_blocks)
    always_tuned(:,iBlock) = all(squeeze(is_tuned(use_tuning_blocks(iBlock),:,:)),2);
end

tuned_cells = all(always_tuned,2);

cell_idx = cell_idx(tuned_cells);

is_sig_diff = zeros(size(tuning_blocks,1),length(cell_idx),length(use_names));
dpd = zeros(size(tuning_blocks,1),length(cell_idx),length(use_names));
for i = 1:length(cell_idx)
    
    unit = cell_idx(i);
    
    for iName = 1:length(use_names)
        n = use_names{iName};
        
        idx_tune_bl = strcmpi({nr(unit).tuning.name},[tw '_' n '_' tuning_blocks{1,1} '_' num2str(tuning_blocks{1,2}(1)) '_' num2str(tuning_blocks{1,2}(2))]);
        
        for iBlock = 1:size(tuning_blocks,1)
            
            % get velocity indices
            idx_tune = strcmpi({nr(unit).tuning.name},[tw '_' n '_' tuning_blocks{iBlock,1} '_' num2str(tuning_blocks{iBlock,2}(1)) '_' num2str(tuning_blocks{iBlock,2}(2))]);
            d = sort( angleDiff(nr(unit).tuning(idx_tune_bl).bootstrap(3,:),nr(unit).tuning(idx_tune).bootstrap(3,:),true,true) );
            % is 0 in the confidence bound? if not, it is different
            if isempty( range_intersection([0 0], [d(ceil(num_iterations*( (1 - conf_level)/2 ))), d(floor(num_iterations*( conf_level + (1-conf_level)/2 )))]) )
                is_sig_diff(iBlock,i,iName) = 1;
            end
            
            idx = find(get_trial_data_indices(trial_data,'date',nr(unit).date));
            if length(trial_data(idx(1)).perturbation_info) > 2
                dir_mult = -1;
            else
                dir_mult = sign(trial_data(idx(1)).perturbation_info(2));
            end
            
            dpd(iBlock,i,iName) = dir_mult*(180/pi)*angleDiff(mean(nr(unit).tuning(idx_tune_bl).bootstrap(3,:),2),mean(nr(unit).tuning(idx_tune).bootstrap(3,:),2),true,true);
        end
    end
end
    
figure;
hold all;
for iName = 1:length(use_names)
    plot(-1,0,'o','Color',plot_colors{iName},'LineWidth',2);
end
legend(use_names);
for iName = 1:length(use_names)
    for iBlock = 1:size(tuning_blocks,1)
        m = mean(squeeze(dpd(iBlock,:,iName)));
        s = 2*std(squeeze(dpd(iBlock,:,iName)))./sqrt(sum(tuned_cells));
        plot(iBlock,m,'o','LineWidth',2,'Color',plot_colors{iName});
        plot([iBlock,iBlock],[m+s, m-s],'LineWidth',2,'Color',plot_colors{iName});
    end
end
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 size(tuning_blocks,1)+1]);

% figure; hist(dpd(logical(all(is_tuned,2)),1),x_bins);
% figure; hist(dpd(logical(all(is_tuned,2)),2),x_bins);

%% Look at data
doAbs = false;

figure;
hold all;
utheta = unique(target_directions);
for i = 1:length(utheta)
    idx = target_directions == utheta(i);
    plot(angleDiff(utheta(i),theta_v(idx),true,~doAbs),'k');
    plot(angleDiff(utheta(i),theta_p(idx),true,~doAbs),'b');
    plot(angleDiff(utheta(i),theta_f(idx),true,~doAbs),'r');
end


figure;
hold all;
utheta = unique(target_directions);
for i = 1:length(utheta)
    idx = target_directions == utheta(i);
    plot(angleDiff(utheta(i),theta_v(idx),true,~doAbs),'b');
    plot(angleDiff(utheta(i),theta_f(idx),true,~doAbs),'r');
    plot(a ngleDiff(theta_f(idx),theta_v(idx),true,~doAbs),'k');
    
end

figure;
hold all;
utheta = unique(target_directions);
for i = 1:length(utheta)
    idx = target_directions == utheta(i);
    plot(angleDiff(utheta(i),theta_f(idx),true,~doAbs),'k');
    plot(angleDiff(utheta(i),theta_v(idx),true,~doAbs),'b');
    plot(angleDiff(utheta(i),theta_f(idx),true,~doAbs),'r');
    
end

%% Histograms of PD changes
d_v = angleDiff(mean(pds_v{1},2),mean(pds_v{2},2),true,true);
d_f = angleDiff(mean(pds_f{1},2),mean(pds_f{2},2),true,true);

close all;
figure; hist(d_f(mean(rs_f{1},2) > 0.1 & mean(rs_v{1},2) > 0.1),20)
figure; hist(d_v(mean(rs_v{1},2) > 0.1 & mean(rs_v{1},2) > 0.1),20)

