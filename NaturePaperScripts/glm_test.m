%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run GLM encoding model
clear; clc; close all;
dataSummary;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session parameters
monkeys = {'Chewie'};
tasks = {'CO'};
perts = {'FF'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trial parameters
num_bins  = 1;
badtrial_params = struct();% ...
%     'ranges', {{'idx_go_cue','idx_movement_on',[5 50]}});
idx_start = {'idx_target_on', 0};
idx_end   = {'idx_trial_end', 0};
train_info = {'AD',[0.5 1]};
test_info = {'AD',[1 80]};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kinematics parameters
kin_rcb_params = struct( ...
    'which_vars', {{'pos','vel','speed'}}, ...
    'rcb_n',      3, ...
    'rcb_hpeaks', [0.01,0.2], ...
    'rcb_b',      0.2, ...
    'flip_time',  true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% array/neuron parameters
in_array  = 'PMd';
out_array = 'M1';
min_fr    = 3; % in Hz
in_dims   = 16;
out_dims  = 8;
spike_rcb_params = struct( ...
    'which_vars', [in_array '_spikes'], ...
    'rcb_n',      0, ...
    'rcb_hpeaks', [0.01,0.2], ...
    'rcb_b',      0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariates parameters
cov_in.basic = { ...
    'pos','all'; ...
    'vel','all'; ...
    'speed','all'; ...
    'pos_rcb','all'; ...
    'vel_rcb','all'; ...
    'speed_rcb','all'};
cov_in.pmd  = cat(1,{[in_array '_spikes'],'all'},cov_in.basic);
cov_in.potent  = cat(1,{[in_array out_array '_potent'],'all'},cov_in.basic);
cov_in.null  = cat(1,{[in_array out_array '_null'],'all'},cov_in.basic);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get files to use
session_idx = find( ...
    ismember(filedb.Task,tasks) & ...
    ismember(filedb.Perturbation,perts) & ...
    ~cellfun(@isempty,filedb.FileNames) & ...
    ismember(filedb.Monkey,monkeys) & ...
    ~(ismember(filedb.Monkey,'Mihili') & datenum(filedb.Date) > datenum('2015-01-01')) & ...
    cellfun(@(x) all(ismember({'M1','PMd'},x)),filedb.Arrays));

[glm_info,pr2,rpr2] = deal(repmat(struct(),1,length(session_idx)));
for iFile = 1:length(session_idx)
    file = session_idx(iFile);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load data
    % load('/Users/mattperich/Data/han_td.mat','trial_data');
    disp(['File ' num2str(iFile) ' of ' num2str(length(session_idx))]);
    fname = [filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat'];
    load(fullfile(rootDir,TDDir,fname),'trial_data');
    [~,td] = getTDidx(trial_data,'result','R');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove bad trials
    td = removeBadTrials(td,badtrial_params);
    td = getSpeed(td);
    % Process kinematics with RCB convolution
    td = convBasisFunc(td,kin_rcb_params);
    % Process spikes with RCB convolution to get history
    td = convBasisFunc(td,spike_rcb_params);
    % trim down trials
    td = truncateAndBin(td,idx_start,idx_end);
    % remove low-firing cells
    td = removeBadNeurons(td,struct('min_fr',min_fr,'do_shunt_check',true));
    
    train_trials = getTDidx(td,'epoch',train_info{1},'range',train_info{2});
    test_trials = getTDidx(td,'epoch',test_info{1},'range',test_info{2});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ESTIMATE DIMENSIONALITY
    if isempty([in_dims out_dims])
        td_temp = appendTDs( ...
            truncateAndBin(td,{'idx_target_on',0},{'idx_target_on',50}), ...
            truncateAndBin(td,{'idx_go_cue',0},{'idx_go_cue',50}));
        in_dims = estimateDimensionality(td_temp,struct('signal',[in_array '_spikes'],'do_smoothing',true));
        out_dims = estimateDimensionality(td_temp,struct('signal',[out_array '_spikes'],'do_smoothing',true));
        clear td_temp
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Re-bin
    td = truncateAndBin(td,num_bins);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get PC projections, potent, and null spaces
    potentNull_params = struct( ...
        'in_signals',[in_array '_spikes'], ...
        'out_signals',[out_array '_spikes'], ...
        'in_dims',in_dims, ...
        'out_dims',out_dims, ...
        'use_trials',getTDidx(td,'epoch',{'BL'}), ...
        'do_smoothing',true);
    [td,potentNull_info] = getPotentSpace(td,potentNull_params);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit the models
    model_names = fieldnames(cov_in);
    for m = 1:length(model_names)
        disp(['Calculating ' model_names{m} ' GLM...']); tic;
        [td, glm_info(iFile).(model_names{m})] = getModel(td,struct( ...
            'model_type','glm', ...
            'in_signals',{cov_in.(model_names{m})}, ...
            'out_signals',{{[out_array '_spikes'],'all'}}, ...
            'model_name',model_names{m}, ...
            'train_idx',train_trials)); toc;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate performance of models
    for m = 1:length(model_names)
        disp(['Evaluating ' model_names{m} ' pR2...']); tic;
        pr2(iFile).(model_names{m}) = evalModel(td,struct( ...
            'model_type','glm', ...
            'out_signals',{glm_info(iFile).(model_names{m}).out_signals}, ...
            'model_name',glm_info(iFile).(model_names{m}).model_name, ...
            'eval_metric','pr2', ...
            'trial_idx',1:length(td))); toc;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate relative performance of models
    temp_names = model_names(~strcmpi(model_names,'basic'));
    for m = 1:length(temp_names)
        disp(['Evaluating ' temp_names{m} ' rpR2...']); tic;
        rpr2(iFile).(temp_names{m}) = evalModel(td,struct( ...
            'model_type','glm', ...
            'out_signals',{glm_info(iFile).(temp_names{m}).out_signals}, ...
            'model_name',{{glm_info(iFile).basic.model_name,glm_info(iFile).(temp_names{m}).model_name}}, ...
            'eval_metric','pr2', ...
            'trial_idx',1:length(td))); toc;
    end, clear temp_names i m;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some plotting
save(fullfile(rootDir,resultsDir,'pmdm1_glm.mat'),'glm_info','pr2','rpr2')


