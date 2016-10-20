function [] = plot_gpfa_3d(trial_data,params_struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot GPFA results
%   vis_data is focused on single trial plotting, including GPFA
%   trajectories, whereas this function will be intended for summarizing
%   across trials
%
% Can also add kinematic traces grouped across trials. Can't add any
%   spiking activity right now because grouping across trials doesn't make
%   sense unless we are smoothing it into a continuous firing rate
%
% Step one: separate out trials by common threads
%   1) target direction
%   2) epoch (BL, AD, WO)
%   3) In the future, task, perturbation type, etc
% Step two: plot trajectories for all trials
%   1) Each dimension against time, color coded by condition
%   2) Any 3 dimensions in trajectory plot, color coded by condition
%
% PARAM_STRUCT OPTIONS:
% Plotting parameters
%   trials               : (vector) trial indices to plot. (Default to all of a given condition)
%   signals              : (cell array) fieldnames of continuous signals to plot (Default to {} for none)
%                            pass in names of fields for continuous signals, e.g. 'vel' or 'acc'
%   plot_conditions      : (cell array) what conditions to plot. Each row represents a condition,
%                            and will be given a unique color in the trajectories.
%                            There are four required columns: {ARRAY, TASK, PERTURBATION, EPOCH}
%                            If you provide it with a string array name, assumes all trials are same condition
%   plot_direction_range : (2 element vector) which range of target directions to plot. Each will get its own columnn. (Defaults to [-pi,pi])
%   dims_for_3d          : dimensions to use for 3D trajectory plot (Default to 1:3)
%   plot_dims            : dimensions to plot against time (Default to first 3 dimensions)
%   gpfa_params          : parameters of GPFA fits (REQUIRED)
%                            bin_width : bin width (in msec) of GPFA model
%                            xdim      : number of assumed latent dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(params_struct,'trials'), trials_to_plot = params_struct.trials; else trials_to_plot = 1:length(trial_data); end
if isfield(params_struct,'signals'), plot_signals = params_struct.signals; else plot_signals = {}; end
if isfield(params_struct,'plot_conditions'), plot_conditions = params_struct.plot_conditions; else plot_conditions = {}; end
if isfield(params_struct,'plot_direction_range'), plot_direction_range = params_struct.plot_direction_range; else plot_direction_range = [-pi,pi]; end
if isfield(params_struct,'plot_dims'), plot_dims = params_struct.plot_dims; else plot_dims = 1:8; end
if isfield(params_struct,'dims_for_3d'), dims_for_3d = params_struct.dims_for_3d; else dims_for_3d = 1:3; end
if isfield(params_struct,'gpfa_params'), gpfa_params = params_struct.gpfa_params; else gpfa_params = []; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   These are parameters that are probably the same all the time
% so it's probably not worth making it an input parameter, but they're here
data_bin_size = 10; %bin size of data in msec
pos_offset = [3-1.7, -33+2.2]; % offset to zero position
event_db = {'idx_target_on','tgt'; ... % list of possible field names for events and a shorthand name
    'idx_go_cue','go'; ...         % add any new events here
    'idx_movement_on','mv'; ...
    'idx_peak_speed','pk'; ...
    'idx_reward','rwd'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   These are a lot of parameters for plotting
% Presumably we won't change these but just in case they are easy to find
font_size   = 12;       % default font size
line_width  = 1;        % standard line width
traj3d_cols = 2;        % how many columns for 3d trajectory plot
time_cols   = 3;        % how many columns for time-variable plots
kin_rows    = 3;        % how many rows for kinematic plots
traj_rows   = 4;        % how many rows for time-varying trajectory plots
pos_location = 'right'; % if position plot is on 'left' or 'right'
trial_event_colors = [0    0.4470    0.7410; ... % using default Matlab r2014b color order for trial events
    0.8500    0.3250    0.0980; ...
    0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(plot_conditions)
    error('Need to specify at least an array to plot in plot_conditions.');
elseif isstring(plot_conditions) % this should be specifying an array
    plot_conditions = {plot_conditions};
end

% find what directions are within plot_direction_range
all_directions = unique([trial_data.target_direction]);
plot_directions = all_directions( all_directions >= plot_direction_range(1) & all_directions <= plot_direction_range(2) );

% split out trials based on parameters of plot_conditions and plot_directions
trial_blocks = cell(size(plot_conditions,1),length(plot_directions));
for iCond = 1:length(plot_conditions)
    for iDir = 1:length(plot_directions)
        trial_blocks{iCond,iDir} = find(get_trial_data_indices(trial_data,'target_direction',plot_directions,'task',plot_conditions{iCond,2},'perturbation',plot_conditions{iCond,3},'epoch',plot_conditions{iCond,4}));
    end
end

% get all available epochs
all_epochs = unique({trial_data.epoch});
all_tasks = unique({trial_data.task});
all_perts = unique({trial_data.perturbation});
all_dirs = unique([trial_data.target_direction]);

% allow for a variable number of events named 'idx_EVENT'
fn = fieldnames(trial_data);
events = fn(cellfun(@(x) ~isempty(regexp(x,'idx_','ONCE')),fn));
clear fn;

% find how many rows are needed
num_rows = length(plot_signals)*kin_rows + length(plot_dims)*traj_rows;

num_cols = traj3d_cols + time_cols;
% use this to partition the subplot space
subplot_grid = repmat((0:num_rows-1)'*num_cols,1,num_cols) + repmat(1:num_cols,num_rows,1);

% some variables to position the columns
switch lower(pos_location)
    case 'left'
        traj3d_start = 0;
        time_start = traj3d_cols;
    case 'right'
        time_start = 0;
        traj3d_start = time_cols;
end

% Make new figure
figure('units','normalized','outerposition',[0.1 0 .85 1]);

% loop along conditions
for iCond = 1:size(plot_conditions,1)
    gpfa_array = plot_conditions{iCond,1};
    
    
    % check that data exists and params provided
    if isfield(trial_data,[gpfa_array '_gpfa'])
        if ~isempty(gpfa_params)
            gpfa_bin_w = gpfa_params.([gpfa_array '_gpfa']).bin_width;
            gpfa_x_dim = gpfa_params.([gpfa_array '_gpfa']).xdim;
        else
            error('GPFA parameter struct input not provided. See documentation.');
        end
    else
        error('GPFA data not present in trial_data.');
    end
    
    % loop along desired target directions
    for iDir = 1:length(plot_directions)
        
        temp_trial_data = trial_data(trial_blocks{iCond,iDir});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot GPFA trajectory path
        subplot(num_rows,num_cols, ...
            reshape(subplot_grid(1:num_rows,traj3d_start+1:traj3d_start+traj3d_cols)',1,traj3d_cols*num_rows));
        hold all;
        for tr_idx = 1:num_trials_to_plot % tr_idx is a dummy variable; useful if you're skipping trials
            tr_num = trials_to_plot(tr_idx); % Use tr_num from here down
            plot3(temp_trial_data(tr_num).([gpfa_array '_gpfa'])(dims_for_3d(1),:), ...
                temp_trial_data(tr_num).([gpfa_array '_gpfa'])(dims_for_3d(2),:), ...
                temp_trial_data(tr_num).([gpfa_array '_gpfa'])(dims_for_3d(3),:), ...
                'linewidth',line_width,'color','k');
            
            % now plot event markers
            for iEvent = 1:length(events)
                % scale bin index to fit gpfa bins
                idx = round(temp_trial_data(tr_num).(events{iEvent})*(data_bin_size/gpfa_bin_w));
                
                plot3(temp_trial_data(tr_num).([gpfa_array '_gpfa'])(dims_for_3d(1),idx), ...
                    temp_trial_data(tr_num).([gpfa_array '_gpfa'])(dims_for_3d(2),idx), ...
                    temp_trial_data(tr_num).([gpfa_array '_gpfa'])(dims_for_3d(3),idx), ...
                    'o','linewidth',3,'color',trial_event_colors(iEvent,:));
            end
        end
        set(gca,'XTick',[],'YTick',[],'ZTick',[], ...
            'Box','off','TickDir','out','FontSize',font_size);
        xlabel(['Factor ' num2str(dims_for_3d(1))],'FontSize',font_size);
        ylabel(['Factor ' num2str(dims_for_3d(2))],'FontSize',font_size);
        zlabel(['Factor ' num2str(dims_for_3d(3))],'FontSize',font_size);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop along the specified continuous signals
        % keep count of how many rows have been used
        row_tally = 0;
        for iSignal = 1:length(plot_signals)
            % Plot kinematics
            subplot(num_rows,num_cols, ...
                reshape(subplot_grid(row_tally+1:row_tally+kin_rows,time_start+1:time_start+time_cols)',1,(num_cols-traj3d_cols)*kin_rows ));
            
            if strcmpi(plot_signals{iSignal},'pos')
                offset = pos_offset;
            else
                offset = [0 0];
            end
            
            hold all;
            for tr_idx = 1:num_trials_to_plot % tr_idx is a dummy variable; useful if you're skipping trials
                tr_num = trials_to_plot(tr_idx); % Use tr_num from here down
                plot(trial_data(tr_num).(plot_signals{iSignal})(:,1) - offset(1),'r','LineWidth',line_width) % x
                plot(trial_data(tr_num).(plot_signals{iSignal})(:,2) - offset(2),'b','LineWidth',line_width) % y
                xlim([1 size(trial_data(tr_num).pos,1)]);
                ylabel(plot_signals{iSignal},'FontSize',font_size);
                set(gca,'Box','off','TickDir','out','XTick',[],'FontSize',font_size);
            end
            row_tally = row_tally + kin_rows;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot GPFA dimension traces over time
        for iDim = 1:length(plot_dims)
            if plot_dims(iDim) <= gpfa_x_dim
                subplot(num_rows,num_cols, ...
                    reshape(subplot_grid(row_tally+1:row_tally+traj_rows,time_start+1:time_start+time_cols)',1,(num_cols-traj3d_cols)*traj_rows ));
                hold all;
                for tr_idx = 1:num_trials_to_plot % tr_idx is a dummy variable; useful if you're skipping trials
                    tr_num = trials_to_plot(tr_idx); % Use tr_num from here down
                    
                    plot(trial_data(tr_num).([gpfa_array '_gpfa'])(plot_dims(iDim),:),'k','LineWidth',line_width);
                end
                axis('tight');
                set(gca,'Box','off','TickDir','out','YTick',[],'XTickLabels',[],'FontSize',font_size);
                ylabel([gpfa_array ' ' num2str(plot_dims(iDim))],'FontSize',font_size)
                
                row_tally = row_tally + traj_rows;
            else
                warning(['Requested dimension (' num2str(plot_dims(iDim)) ') is larger than available latent dimensions (xDim = ' num2str(gpfa_x_dim) '. Skipping.']);
            end
        end
    end
    
    xlabel('Time (bins)','FontSize',font_size);
end

end
