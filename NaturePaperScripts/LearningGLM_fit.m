%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run GLM encoding model
clear; clc; close all;
dataSummary;

sessions      = { ...
    'Chewie','2016-09-09'; ... % VR
    'Chewie','2016-09-12'; ...
    'Chewie','2016-09-14'; ...
    'Chewie','2016-10-06'; ...
%     'Mihili','2014-03-03'; ...
%     'Mihili','2014-03-04'; ...
%     'Mihili','2014-03-06'; ...
%     'Chewie','2016-09-15'; ... % CF
%     'Chewie','2016-10-05'; ...
%     'Chewie','2016-10-07'; ...
%     'Chewie','2016-10-11'; ...
%     'Mihili','2014-02-03'; ...
%     'Mihili','2014-02-17'; ...
%     'Mihili','2014-02-18'; ...
%     'Mihili','2014-03-07'; ...
    % 'MrT','2013-08-19'; ... % CF
    % 'MrT','2013-08-21'; ...
    % 'MrT','2013-08-23'; ...
    % 'MrT','2013-09-03'; ... %VR
    % 'MrT','2013-09-05'; ...
    % 'MrT','2013-09-09'; ...
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra info for path name
name_info = 'VR';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session parameters
monkeys      = unique(sessions(:,1));
dates        = sessions(:,2);
tasks        = {'CO'};
perts        = {'FF','VR'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trial parameters
trial_params;
num_bins     = 5;
train_info   = {'AD',[0.5 1]};
num_boots    = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kinematics parameters
kin_lags     = 2;
kin_rcb_params = struct( ...
    'which_vars', {{'vel','speed','acc'}}, ...
    'rcb_n',      2, ...
    'rcb_hpeaks', [0.01*num_bins,0.1], ...
    'rcb_b',      0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% array/neuron parameters
in_array     = 'PMd';
out_array    = 'M1';
in_dims      = 16;
out_dims     = 8;
null_size    = ''; %'double',''
spike_rcb_params = struct( ...
    'which_vars', [in_array '_spikes'], ...
    'rcb_n',      0, ...
    'rcb_hpeaks', [0.01*num_bins,0.1], ...
    'rcb_b',      0.3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM parameters
do_lasso     = false;
lasso_lambda = 0.05;
lasso_alpha  = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariates parameters
basic_inputs = {'pos','vel','speed','pos_rcb','vel_rcb','speed_rcb'};
basic_inputs = [basic_inputs', repmat({'all'},length(basic_inputs),1)];
cov_in.basic = basic_inputs;

% basic_inputs = {};
cov_in.(lower(in_array))  = cat(1,{[in_array '_spikes'],'all'},basic_inputs);
if ~strcmp(in_array,out_array) % potent/null doesn't make sense if they're the same
%     cov_in.([lower(in_array) '_pca'])  = cat(1,{[in_array '_pca'],1:in_dims},basic_inputs);
%     cov_in.potent  = cat(1,{[in_array out_array '_potent'],'all'},basic_inputs);
%     cov_in.null  = cat(1,{[in_array out_array '_null'],'all'},basic_inputs);
end

ANALYSIS_NAME = [in_array out_array '_glm_' name_info];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build function call list for loadTDfiles
func_calls = [trial_func_calls, { ...
    {@removeBadNeurons,badneuron_params}, ...
    {@getSpeed}, ...
    {@binTD,num_bins}, ...
    {@convBasisFunc,kin_rcb_params}, ...
    {@convBasisFunc,spike_rcb_params}}];

session_idx = getFileDBidx(filedb, ...
    {'Task',tasks,'Perturbation',perts,'Monkey',monkeys,'Date',dates}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});


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
    
    if any(ismember(model_names,{'potent','null','m1_pca','pmd_pca'}))
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
            trial_data(trial).PMd_spikes = td(trial).PMd_spikes;
        end, clear td;
    end
    
    %%% Do that thing where I lag kinematic signals
    if kin_lags > 0
    for trial = 1:length(trial_data)
        for i = 1:length(basic_inputs)
            temp = trial_data(trial).(basic_inputs{i});
            trial_data(trial).(basic_inputs{i}) = [temp(kin_lags+1:end,:); NaN(kin_lags,size(temp,2))];
        end
    end
    end
    
    trial_data = trimTD(trial_data,idx_start,idx_end);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For now we only study AD trials
%     if strcmpi(train_info{1},'ad')
%         disp('Cutting down to AD trials only!');
%         [~,trial_data] = getTDidx(trial_data,'epoch','AD');
%         train_trials = getTDidx(trial_data,'epoch','AD','range',train_info{2});
%     end
    
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
    % Cross validate the models
    
    
    params.glm_info = glm_info;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save it
    save(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[filenames{iFile} '_fit.mat']),'trial_data','params')
end

