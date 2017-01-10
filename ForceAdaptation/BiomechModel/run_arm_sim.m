close all;
clear;
clc;

filenames = { ...
    'Chewie_CO_FF_2013-10-22', ...
    'Chewie_CO_FF_2013-10-23', ...
    'Chewie_CO_FF_2013-10-31', ...
    'Chewie_CO_FF_2013-11-01', ...
    'Chewie_CO_FF_2013-12-03', ...
    'Chewie_CO_FF_2013-12-04', ...
    'Chewie_CO_FF_2015-07-01', ...
    'Chewie_CO_FF_2015-07-03', ...
    'Mihili_CO_FF_2014-02-03', ...
    'Mihili_CO_FF_2014-02-17', ...
    'Mihili_CO_FF_2014-02-18', ...
    'Mihili_CO_FF_2014-03-07', ...
    'Mihili_CO_FF_2015-06-11', ...
    'Mihili_CO_FF_2015-06-17', ...
    'Chewie_CO_FF_2015-06-29', ...
    'Chewie_CO_FF_2015-06-30', ...
    'Chewie_CO_FF_2015-07-06', ...
    'Chewie_CO_FF_2015-07-07', ...
    'Chewie_CO_FF_2015-07-08', ...
    'Mihili_CO_FF_2015-06-10', ...
    'Mihili_CO_FF_2015-06-15', ...
    'Mihili_CO_FF_2015-06-16', ...
    };

for iFile = 1:length(filenames)
    
    clearvars -except filenames iFile;
    close all;
    clc;
    iFile
    
    root_dir = 'F:\TrialDataFiles\';
    out_dir = fullfile(root_dir,'biomech_model\');
    do_sim = false; %if false just redo tuning and use old model
    do_neurons = true;
    do_tuning = true;
    
    do_sliding_window_bullshit = true;
    tune_real_neurons = false;
    
    % load in trial data
    load([root_dir filenames{iFile} '.mat'],'trial_data');
    
    %
    use_models = {'muscle'};
    left_arm = true; % flips from right arm to left arm
    
    %%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%
    pos_offset = [1,-31]; % position offset from behavior
    origin_pos = [15 25];
    dt = 0.01;
    
    %%% biomech model
    % use Cheng/Scott 2000 for parameters
    switch lower(trial_data(1).monkey)
        case 'chewie'
            body_weight = 12; % body mass in kg 11.7
            M1 = 34.4*body_weight/1000; % from Scott 2000
            M2 = 25.2*body_weight/1000;
            L1 = 0.25;
            L2 = 0.22;
            %             M1 = 0.7; % mass of upper arm in kg
            %             M2 = 0.7; % mass of lower arm in kg
            %             L1 = 0.17; % length of upper arm in m
            %             L2 = 0.17; % length of lower arm in m
        case 'mihili'
            body_weight = 10; %8.8
            M1 = 34.4*body_weight/1000;
            M2 = 25.2*body_weight/1000;
            L1 = 0.20;
            L2 = 0.24;
            %             M1 = 0.55; % mass of upper arm in kg
            %             M2 = 0.55; % mass of lower arm in kg
            %             L1 = 0.14; % length of upper arm in m
            %             L2 = 0.14; % length of lower arm in m
    end
    
    TH_1_min = 0;  % minimum shoulder angle
    TH_1_max = pi; % maximum shoulder angle
    TH_2_min = 0;  % minimum elbow angle
    TH_2_max = pi; % maximum elbow angle
    
    %%% muscle model
    simple_muscle_model = false;
    use_insertion_angles = false;
    
    %%% neural activity model
    num_neurons = 500;
    muscle_gains = [1,1,1,1,1,1]; % gain terms for [sh flex, sh ext, el flex, el ext]
    use_muscle_model = 'all'; %'single','joint','synergy','all'
    weight_distribution = 'gaussian'; %'uniform','gaussian','skewgauss'
    mean_lag = 10; % mean lag in bins, that spikes will be shifted before tuning
    std_lag = 5; % how many bins of standard deviation
    
    if tune_real_neurons
        mean_lag = 10;
    end
    
    if ~iscell(use_models)
        use_models = {use_models};
    end
    
    % build params struct
    params = struct('use_models',use_models,'pos_offset',pos_offset,'origin_pos',origin_pos,'dt',dt,'M1',M1,'M2',M2,'L1',L1,'L2',L2, ...
        'num_neurons',num_neurons,'muscle_gains',muscle_gains,'use_muscle_model',use_muscle_model,'mean_lag',mean_lag,'std_lag',std_lag, ...
        'simple_muscle_model',simple_muscle_model,'weight_distribution',weight_distribution,'left_arm',left_arm);
    
    if isfield(trial_data(1),'result'),
        trial_data = trial_data(getTDidx(trial_data,'result','R'));
    end
    
    %trial_data = trial_data(getTDidx(trial_data,'epoch','BL'));
        
    if do_sim
        % do modeling
        disp('Calculating kinematics and dynamics...');
        
        arm_sim_muscles;

        % save data
        save([out_dir filenames{iFile} '_sim.mat'],'sim_data','params');
    end
    
    if do_neurons
        % Now, generate neural activity for each trial
        % save data
        load([out_dir filenames{iFile} '_sim.mat'],'sim_data','params');
        
        arm_sim_neurons;

        % save data
        save([out_dir filenames{iFile} '_neurons.mat'],'sim_data','neural_tcs','params');
        
    end
    
    if do_tuning
        disp('Using old data...')
        load([out_dir filenames{iFile} '_neurons.mat'],'sim_data','neural_tcs','params');
        
        arm_sim_tuning;
        
        save([out_dir filenames{iFile} '_tuning.mat'],'sw_data','tc_data','params');
        
    end
end

%% Plot the simulation
if 1
    %     clear params;
    close all;
    %     params.dt = dt;
    %     params.M1 = M1; % mass of upper arm in kg
    %     params.M2 = M2; % mass of lower arm in kg
    %     params.L1 = L1; % length of upper arm in m
    %     params.L2 = L2; % length of lower arm in m
    %     params.pos_offset = pos_offset; % position offset from behavior
    %     params.origin_pos = origin_pos;
    %     params.dt = dt;
    %     params.M_max = M_max;
    %     params.M_min = M_min;
    %     params.T_max = T_max;
    %     params.T_min = T_min;
    %     params.num_neurons = num_neurons;
    %     params.comp_blocks = [1 2 3];
    %
    params.type = 'time_signals';
%             params.type = 'freeze_video';
    params.resolution = 5;
    params.signals = {'vel','torques','muscles','muscle_neurons'};
    params.kin_model = 'real';
    params.comp_blocks = [1,2,3];
    
    %     params.torques_lim = [-0.5,0.25];
    %     params.vel_lim = [-0.3,0.1];
    
    I=cell2mat({trial_data.target_direction})==trial_data(58).target_direction & strcmpi({trial_data.epoch},'AD'); idx = find(I);
    %         idx = 58; % baseline FOR MIHILI 2015-06-16
    %     idx = 348; % decent ad in same direction

%     idx = [10,211];
% idx = find(getTDidx(trial_data,'epoch','BL'));
%     [~,I] = sort(cellfun(@(x) max(max(x)) - min(min(x)),{sim_data(idx).torques},'Uni',1),2,'ascend');
%     idx = idx(I);

idx = [58,348];
    
    for iTrial = idx
        iTrial
        params.idx_start = trial_data(iTrial).idx_go_cue;
        params.idx_end = trial_data(iTrial).idx_trial_end;
        plot_arm_sim(sim_data(iTrial),params,tc_data);
        
        if length(idx) > 4
            pause;
            close all;
        end
    end
end
