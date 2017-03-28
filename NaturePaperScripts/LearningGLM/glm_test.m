 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run GLM encoding model
clear; clc; close all;
dataSummary;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name for saving
ANALYSIS_NAME = 'PMdM1_glm_machens_50ms';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session parameters
monkeys = {'MrT'};
tasks = {'CO'};
perts = {'FF'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trial parameters
trial_params;
num_bins  = 5;
train_info = {'AD',[0.5 1]};
num_boots = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kinematics parameters
kin_rcb_params = struct( ...
    'which_vars', {{'pos','vel'}}, ...
    'rcb_n',      2, ...
    'rcb_hpeaks', [0.01*num_bins,0.1], ...
    'rcb_b',      0.1, ...
    'flip_time',  true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% array/neuron parameters
in_array  = 'PMd';
out_array = 'M1';
in_dims   = [];
out_dims  = [];
null_size = 'double'; %'double','real'
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
    'vel','all'; ...
    'vel_rcb','all'};
cov_in.pmd  = cat(1,{[in_array '_spikes'],'all'},cov_in.basic);
% cov_in.pmd_pca  = cat(1,{[in_array '_pca'],1:15},cov_in.basic);
% cov_in.potent  = cat(1,{[in_array out_array '_potent'],'all'},cov_in.basic);
% cov_in.null  = cat(1,{[in_array out_array '_null'],'all'},cov_in.basic);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build function call list for loadTDfiles
func_calls = [trial_func_calls, { ...
    {@removeBadNeurons,badneuron_params}, ...
    {@getSpeed}, ...
    {@convBasisFunc,kin_rcb_params}, ...
    {@convBasisFunc,spike_rcb_params}, ...
    {@trimTD,idx_start,idx_end}}];

session_idx = getFileDBidx(filedb, ...
    {'Task',tasks,'Perturbation',perts,'Monkey',monkeys}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
filenames = cell(1,length(session_idx));
for iFile = 1:length(session_idx)
    file = session_idx(iFile);
    filenames{iFile} = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
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
    [trial_data,params] = loadTDfiles(filenames{iFile},func_calls{:});
    train_trials = getTDidx(trial_data,'epoch',train_info{1},'range',train_info{2});
    
    if size(trial_data(1).([in_array '_spikes']),2) == 0 || size(trial_data(1).([out_array '_spikes']),2) == 0
        error('No signals to predict or use as inputs!');
    end
    
    % get dimensionality if desired
    if isempty([in_dims, out_dims])
        trim_idx = {'idx_movement_on',-50;'idx_movement_on',50};
        [~,td_temp] = getTDidx(trial_data,'epoch','BL');
        td_temp = smoothSignals(td_temp,struct('signals',{{'M1_spikes','PMd_spikes'}},'kernel_SD',pn_kernel_SD));
        temp_in_dims = estimateDimensionality(td_temp,struct( ...
            'signals',[in_array '_spikes'], ...
            'trim_idx',{trim_idx}));
        temp_out_dims = estimateDimensionality(td_temp,struct( ...
            'signals',[out_array '_spikes'], ...
            'trim_idx',{trim_idx}));
        
        switch lower(null_size)
            case 'double'
                in_dims = 2*out_dims;
        end
    else
        temp_in_dims = in_dims;
        temp_out_dims = out_dims;
    end
    
    potentNull_params = struct( ...
        'in_signals',[in_array '_spikes'], ...
        'out_signals',[out_array '_spikes'], ...
        'in_dims',temp_in_dims, ...
        'out_dims',temp_out_dims, ...
        'use_trials',{{'epoch','BL'}});
    trial_data = getPotentSpace(trial_data,potentNull_params);
    trial_data = binTD(trial_data,num_bins);
    params.potentNull_params = potentNull_params;

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate performance of models
    for m = 1:length(model_names)
        disp(['Evaluating ' model_names{m} ' pR2...']); tic;
        pr2 = evalModel(trial_data,struct( ...
            'model_type','glm', ...
            'out_signals',{glm_info.(model_names{m}).out_signals}, ...
            'model_name',glm_info.(model_names{m}).model_name, ...
            'eval_metric','pr2', ...
            'trial_idx',1:length(trial_data), ...
            'num_boots',num_boots, ...
            'do_parallel',true)); toc;
                % add it to the output struct
        for trial = 1:size(pr2,1)
            trial_data(trial).([model_names{m} '_pr2']) = squeeze(pr2(trial,:,:));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate relative performance of models
    temp_names = model_names(~strcmpi(model_names,'basic'));
    for m = 1:length(temp_names)
        disp(['Evaluating ' temp_names{m} ' rpR2...']); tic;
        rpr2 = evalModel(trial_data,struct( ...
            'model_type','glm', ...
            'out_signals',{glm_info.(temp_names{m}).out_signals}, ...
            'model_name',{{glm_info.basic.model_name,glm_info.(temp_names{m}).model_name}}, ...
            'eval_metric','pr2', ...
            'trial_idx',1:length(trial_data), ...
            'num_boots',num_boots, ...
            'do_parallel',true)); toc;
        % add it to the output struct
        for trial = 1:size(rpr2,1)
            trial_data(trial).([temp_names{m} '_rpr2']) = squeeze(rpr2(trial,:,:));
        end
    end, clear temp_names i m;
    
    params.glm_info = glm_info;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save it
    [~,fname,~] = fileparts(filenames{iFile});
    save(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[fname '.mat']),'trial_data','params')
end


%%
clc; close all; clearvars -except ANALYSIS_NAME;
dataSummary;
load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,'start_params.mat'));

models = {'pmd'};%,'potent','null'};

figure; hold all;
for m = 1:length(models)
    a=[];
    for iFile = 1:length(filenames)
        [~,fname,~] = fileparts(filenames{iFile});
        load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[fname '.mat']));
    
        train_idx = params.glm_info.(models{m}).train_idx;
        ad_idx = getTDidx(trial_data,'epoch','AD','range',[1 100]);
        
        temp = squeeze(mean(cat(3,trial_data.([models{m} '_rpr2'])),2))';
        
        cv_ref = mean(temp(train_idx(1)-15:train_idx(1)-1,:),1);
        idx = cv_ref > 0.01;
        
        temp = moving_average(temp,10);
        
        temp = (temp-repmat(cv_ref,size(temp,1),1))./repmat(cv_ref,size(temp,1),1);
        a = [a, temp(ad_idx,idx)];
%         a = [a, temp(train_idx(1)-150:train_idx(1)+50,idx)];
        
    end
    
    temp = median(a,2);
    plot(temp);
    % plot([train_idx(1) train_idx(1)]-1,[-1 1],'k--');
    % plot([train_idx(end) train_idx(end)]-1,[-1 1],'k--');
    
end


