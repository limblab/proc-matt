%% Combine data files and build new structure
dataSummary;

idx = strcmpi(sessionList(:,1),'Chewie');% & strcmpi(sessionList(:,4),'CO');
sessionList = sessionList(idx,:);

for iFile = 1:size(sessionList,1)
    disp(['File ' num2str(iFile) ' of ' num2str(size(sessionList,1))]);
    
    data = combineFiles({['F:\' sessionList{iFile,1} '\Processed\' sessionList{iFile,2} '\' sessionList{iFile,4} '_' sessionList{iFile,3} '_BL_' sessionList{iFile,2} '.mat'], ...
        ['F:\' sessionList{iFile,1} '\Processed\' sessionList{iFile,2} '\' sessionList{iFile,4} '_' sessionList{iFile,3} '_AD_' sessionList{iFile,2} '.mat'], ...
        ['F:\' sessionList{iFile,1} '\Processed\' sessionList{iFile,2} '\' sessionList{iFile,4} '_' sessionList{iFile,3} '_WO_' sessionList{iFile,2} '.mat']});
    
    trial_data = parseFileByTrial(data,struct('dt',0.0005));
    
    save(['F:\trial_data_files\' sessionList{iFile,1} '_' sessionList{iFile,4} '_' sessionList{iFile,3} '_' sessionList{iFile,2} '.mat'],'trial_data');
end
disp('Done.');

%%
% % % % Load in data and build one giant master array
% % % dataSummary;
% % % all_trial_data = [];
% % % for iFile = 1:size(sessionList,1)
% % %     load(['F:\trial_data_files\' sessionList{iFile,1} '_' sessionList{iFile,4} '_' sessionList{iFile,3} '_' sessionList{iFile,2} '.mat'],'trial_data');
% % %     
% % %     if ~isfield(trial_data,'M1_spikes');
% % %         [trial_data.M1_spikes] = deal([]);
% % %     end
% % %     if ~isfield(trial_data,'PMd_spikes');
% % %         [trial_data.PMd_spikes] = deal([]);
% % %     end
% % %     
% % %     all_trial_data = [all_trial_data trial_data];
% % % end
% % % 
% % % clear trial_data;
% % % trial_data = all_trial_data;
% % % save('F:\trial_data_files\All_CO_FF.mat','trial_data')
% % % 
% % % %% Test building GPFA model
% % % root_dir = 'F:\trial_data_files\';
% % % % filename = 'Chewie_CO_FF_2015-07-07';
% % % filename = 'Mihili_CO_FF_2014-02-17';
% % % load(fullfile(root_dir,[filename '.mat']));
% % % 
% % % % try the three different variants to make sure they work
% % % params.bin_w = 30;
% % % params.xdim = 8;
% % % params.kernsd = 30;
% % % params.arrays = 'M1';
% % % [trial_data, M1_gpfa] = run_gpfa(trial_data,params);
% % % 
% % % % params.arrays = 'PMd';
% % % % [trial_data, PMd_gpfa] = run_gpfa(trial_data,params);
% % % % 
% % % % params.arrays = 'Both';
% % % % [trial_data, Both_gpfa] = run_gpfa(trial_data,params);
% % % 
% % % % store the models for later use
% % % gpfa_models.M1_gpfa = M1_gpfa.model;
% % % % gpfa_models.PMd_gpfa = PMd_gpfa.model;
% % % % gpfa_models.Both_gpfa = Both_gpfa.model;
% % % gpfa_params.M1_gpfa = M1_gpfa.params;
% % % % gpfa_params.PMd_gpfa = PMd_gpfa.params;
% % % % gpfa_params.Both_gpfa = Both_gpfa.params;
% % % 
% % % % save in a new file
% % % save(fullfile(root_dir,[filename '_gpfa.mat']),'trial_data','gpfa_models','gpfa_params');
% % % 
% % % %% Test vis_data
% % % % close all;
% % % clear params;
% % % params.plot_gpfa = false;
% % % params.gpfa_dims = 1:3;
% % % params.gpfa_array = 'M1';
% % % params.trials = 190;
% % % params.signals = {'vel','acc'};
% % % % params.gpfa_params = gpfa_params;
% % % vis_data(trial_data,params);
% % % 
% % % %% Test plot_gpfa
% % % load('F:\trial_data_files\Mihili_CO_FF_2014-02-17_gpfa.mat');
% % % % load('F:\trial_data_files\Chewie_CO_FF_2015-07-07_gpfa.mat');
% % % % trim trials with empty peak speed idx
% % % trial_data = trial_data(~cellfun(@isempty,{trial_data.idx_peak_speed}));
% % % 
% % % % Prune bad trials
% % % idx = find(get_trial_data_indices(trial_data,'epoch','BL'));
% % % bl_reaction_times = [trial_data(idx).idx_movement_on]-[trial_data(idx).idx_go_cue];
% % % bl_peak_times = [trial_data(idx).idx_peak_speed]-[trial_data(idx).idx_movement_on];
% % % bl_trial_times = [trial_data(idx).idx_reward]-[trial_data(idx).idx_movement_on];
% % % 
% % % idx = find(get_trial_data_indices(trial_data,'epoch','AD'));
% % % ad_reaction_times = [trial_data(idx).idx_movement_on]-[trial_data(idx).idx_go_cue];
% % % ad_peak_times = [trial_data(idx).idx_peak_speed]-[trial_data(idx).idx_movement_on];
% % % ad_trial_times = [trial_data(idx).idx_reward]-[trial_data(idx).idx_movement_on];
% % % 
% % % bl_bad_trials = find( (bl_reaction_times <= mean(bl_reaction_times) - 3*std(bl_reaction_times)) | ...
% % %     (bl_reaction_times >= mean(bl_reaction_times) + 3*std(bl_reaction_times)) | ...
% % %     (bl_peak_times <= mean(bl_peak_times) - 3*std(bl_peak_times)) | ...
% % %     (bl_peak_times >= mean(bl_peak_times) + 3*std(bl_peak_times)) | ...
% % %     (bl_trial_times <= mean(bl_trial_times) - 3*std(bl_trial_times)) | ...
% % %     (bl_trial_times >= mean(bl_trial_times) + 3*std(bl_trial_times)) );
% % % ad_bad_trials = find( (ad_reaction_times <= mean(bl_reaction_times) - 3*std(bl_reaction_times)) | ...
% % %     (ad_reaction_times >= mean(bl_reaction_times) + 3*std(bl_reaction_times)) | ...
% % %     (ad_peak_times <= mean(bl_peak_times) - 3*std(bl_peak_times)) | ...
% % %     (ad_peak_times >= mean(bl_peak_times) + 3*std(bl_peak_times)) | ...
% % %     (ad_trial_times <= mean(bl_trial_times) - 3*std(bl_trial_times)) | ...
% % %     (ad_trial_times >= mean(bl_trial_times) + 3*std(bl_trial_times)) );
% % % 
% % % close all;
% % % clc;
% % % clear params;
% % % params.gpfa_array = 'M1';
% % % params.plot_dims = 1:3;
% % % params.signals = {'vel_x','vel_y'};
% % % idx_bl = find(get_trial_data_indices(trial_data,'epoch','BL'));
% % % idx_bl(bl_bad_trials) = [];
% % % idx_ad = find(get_trial_data_indices(trial_data,'epoch','AD'));
% % % idx_ad(ad_bad_trials) = [];
% % % params.trial_conditions = {idx_bl, idx_ad(1:floor(length(idx_ad)/2)), idx_ad(end-ceil(length(idx_ad)/2):end)};
% % % params.max_trials = 100;
% % % params.gpfa_params = gpfa_params;
% % % params.plot_direction_range = [0 0];
% % % params.plot_3d = true;
% % % 
% % % params.align_idx = {'idx_target_on',0; 'idx_go_cue',0};
% % % plot_gpfa(trial_data,params);
% % % 
% % % params.align_idx = {'idx_go_cue',0; 'idx_peak_speed',0};
% % % plot_gpfa(trial_data,params);
% % % 
% % % params.align_idx = {'idx_peak_speed',0; 'idx_reward',-400};
% % % plot_gpfa(trial_data,params);
% % % 
% % % 
% % % 
