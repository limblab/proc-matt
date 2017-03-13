% DIFFERENCES IN OLD VERSION
% make sure FR is always non-zero
% do FR min only on baseline
% I was square rooting for PCA, but not projecting square roots, I think
% I was smoothing for PCA, but not projecting smoothed (I think?)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run GLM encoding model
clear; clc; close all;
dataSummary;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name for saving
ANALYSIS_NAME = 'PMdM1_glm_test';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session parameters
monkeys = {'Chewie'};
tasks = {'CO'};
perts = {'FF'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trial parameters
trial_params;
num_bins  = 5;
train_info = {'AD',[0.5 1]};
num_boots = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kinematics parameters
kin_rcb_params = struct( ...
    'which_vars', {{'vel','speed'}}, ...
    'rcb_n',      2, ...
    'rcb_hpeaks', [0.01*num_bins,0.1], ...
    'rcb_b',      0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% array/neuron parameters
in_array  = 'PMd';
out_array = 'M1';
in_dims   = 16;
out_dims  = 8;
null_size = 'double'; %'double',''
spike_rcb_params = struct( ...
    'which_vars', [in_array '_spikes'], ...
    'rcb_n',      0, ...
    'rcb_hpeaks', [0.01*num_bins,0.1], ...
    'rcb_b',      0.3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM parameters
do_lasso = false;
lasso_lambda = 0.05;
lasso_alpha = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariates parameters
cov_in.basic = { ...
    'pos'; 'vel'; 'speed'; 'vel_rcb';'speed_rcb'};
cov_in.basic = [cov_in.basic, repmat({'all'},size(cov_in.basic,1),1)];
% cov_in.pmd  = cat(1,{[in_array '_spikes'],'all'},cov_in.basic);
% cov_in.pmd_pca  = cat(1,{[in_array '_pca'],1:15},cov_in.basic);
cov_in.potent  = cat(1,{[in_array out_array '_potent'],'all'},cov_in.basic);
cov_in.null  = cat(1,{[in_array out_array '_null'],'all'},cov_in.basic);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build function call list for loadTDfiles
func_calls = [trial_func_calls, { ...
    {@removeBadNeurons,badneuron_params}, ...
    {@getSpeed}, ...
    {@binTD,num_bins}, ...
    {@convBasisFunc,kin_rcb_params}, ...
    {@convBasisFunc,spike_rcb_params}, ...
    {@trimTD,idx_start,idx_end}}];

session_idx = getFileDBidx(filedb, ...
    {'Task',tasks,'Perturbation',perts,'Monkey',monkeys}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});

session_idx = session_idx(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
[filepaths,filenames] = deal(cell(1,length(session_idx)));
for iFile = 1:length(session_idx)
    file = session_idx(iFile);
    filepaths{iFile} = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
    filenames{iFile} = [filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file}];
end

if ~exist(fullfile(rootDir,resultsDir,ANALYSIS_NAME),'dir')
    mkdir(fullfile(rootDir,resultsDir,ANALYSIS_NAME))
else
    disp('Warning: will overwrite existing results.');
end
save(fullfile(rootDir,resultsDir,ANALYSIS_NAME,'start_params.mat'));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run GLM
model_names = fieldnames(cov_in);
for iFile = 1:length(session_idx)
    disp(['Fitting GLMs for File ' num2str(iFile) ' of ' num2str(length(session_idx)) '.']);
    file = session_idx(iFile);
    
    % get relevant trials
    [trial_data,params] = loadTDfiles(filepaths{iFile},func_calls{:});
    train_trials = getTDidx(trial_data,'epoch',train_info{1},'range',train_info{2});
    
    if size(trial_data(1).([in_array '_spikes']),2) == 0 || size(trial_data(1).([out_array '_spikes']),2) == 0
        error('No signals to predict or use as inputs!');
    end
    
    if length(model_names) > 1
    % get dimensionality if desired
    if isempty([in_dims, out_dims])
        disp('Getting dimensionality with Machens method.');
        switch lower(null_size)
            case 'double'
                [~,td_temp] = getTDidx(trial_data,'epoch',{'BL','AD'});
                td_temp = smoothSignals(td_temp,struct('signals',{{'M1_spikes','PMd_spikes'}},'kernel_SD',pn_kernel_SD));
                td_temp = trimTD(td_temp,{'idx_movement_on',-round(30/num_bins)},{'idx_movement_on',round(50/num_bins)});
                temp_out_dims = estimateDimensionality(td_temp,struct( ...
                    'signals',[out_array '_spikes']));
                temp_in_dims = 2*temp_out_dims;
            otherwise
                %trim_idx = {'idx_movement_on',-50;'idx_movement_on',50};
                [~,td_temp] = getTDidx(trial_data,'epoch',{'BL','AD'});
                td_temp = smoothSignals(td_temp,struct('signals',{{'M1_spikes','PMd_spikes'}},'kernel_SD',pn_kernel_SD));
                td_temp = appendTDs( ...
                    trimTD(td_temp,{'idx_target_on',0},{'idx_target_on',round(50/num_bins)}), ...
                    trimTD(td_temp,{'idx_go_cue',0},{'idx_go_cue',round(15/num_bins)}), ...
                    trimTD(td_temp,{'idx_movement_on',0},{'idx_movement_on',round(50/num_bins)}));
                temp_in_dims = estimateDimensionality(td_temp,struct( ...
                    'signals',[in_array '_spikes']));
                temp_out_dims = estimateDimensionality(td_temp,struct( ...
                    'signals',[out_array '_spikes']));
        end, clear td_temp;
    else
        temp_in_dims = in_dims;
        temp_out_dims = out_dims;
    end
    
    disp('Doing some smoothing shit')
    
    td = trial_data;
    
        trial_data = sqrtTransform(trial_data,{'M1_spikes','PMd_spikes'});
    trial_data = smoothSignals(trial_data,struct('signals',{{'M1_spikes','PMd_spikes'}},'kernel_SD',0.1));
    
    potentNull_params = struct( ...
        'in_signals',[in_array '_spikes'], ...
        'out_signals',[out_array '_spikes'], ...
        'in_dims',temp_in_dims, ...
        'out_dims',temp_out_dims, ...
        'use_trials',{{'epoch','BL'}}, ...
        'sqrt_transform',false, ...
        'do_smoothing',false, ...
        'kernel_SD',pn_kernel_SD);
    trial_data = getPotentSpace(trial_data,potentNull_params);
    params.potentNull_params = potentNull_params;
    
    for trial = 1:length(td)
        trial_data(trial).M1_spikes = td(trial).M1_spikes;
    end, clear td;

    end
    
%     trial_data = fuck_this(trial_data);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For now we only study AD trials
    if 1
    disp('Cutting down to AD trials only!');
    [~,trial_data] = getTDidx(trial_data,'epoch','AD');
    train_trials = getTDidx(trial_data,'epoch','AD','range',train_info{2});
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit the models
    for m = 1:length(model_names)
        disp(['Calculating ' model_names{m} ' GLM...']); tic;
        [trial_data, glm_info.(model_names{m})] = getModel(trial_data,struct( ...
            'model_type','glm', ...
            'in_signals',{cov_in.(model_names{m})}, ...
            'out_signals',{{[out_array '_spikes'],'all'}}, ...
            'model_name',model_names{m}, ...
            'train_idx',train_trials, ...
            'do_lasso',do_lasso, ...
            'lasso_lambda',lasso_lambda, ...
            'lasso_alpha',lasso_alpha)); toc;
    end
    
    params.glm_info = glm_info;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save it
    save(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[filenames{iFile} '_fit.mat']),'trial_data','params')
end

%%
LearningGLM_eval;

