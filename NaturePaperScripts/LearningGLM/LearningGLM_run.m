%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run GLM encoding model
clear; clc; close all;
dataSummary;

sessions      = { ...
    'Chewie','2016-09-09'; ... % VR
    'Chewie','2016-09-12'; ...
    'Chewie','2016-09-14'; ...
    'Chewie','2016-10-06'; ...
    'Mihili','2014-03-03'; ...
    'Mihili','2014-03-04'; ...
    'Mihili','2014-03-06'; ...
    'Chewie','2016-09-15'; ... % CF
    'Chewie','2016-10-05'; ...
    'Chewie','2016-10-07'; ...
    'Chewie','2016-10-11'; ...
    'Mihili','2014-02-03'; ...
    'Mihili','2014-02-17'; ...
    'Mihili','2014-02-18'; ...
    'Mihili','2014-03-07'; ...
%     'MrT','2013-08-19'; ... % CF
%     'MrT','2013-08-21'; ...
%     'MrT','2013-08-23'; ...
%     'MrT','2013-09-03'; ... %VR
%     'MrT','2013-09-05'; ...
%     'MrT','2013-09-09'; ...
    };


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra info for path name
name_info = 'final';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session parameters
monkeys      = unique(sessions(:,1));
dates        = sessions(:,2);
tasks        = {'CO'};
perts        = {'FF'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trial parameters
num_bins     = 5;
use_results = {'R'};
idx_start = {'idx_trial_start', 0};
idx_end   = {'idx_trial_end', -2};
badtrial_params = struct();%'ranges',{{'idx_go_cue','idx_movement_on',[5 50]}});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM parameters
train_info   = {'AD',[0.5 1]};
cv_blocks    = 10; % if > 0, cross validates on training data
cv_block_size = 10;
cv_type      = 'fold'; %'rand','fold'
num_boots    = 0;
do_lasso     = false;
lasso_lambda = 0.05;
lasso_alpha  = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kinematics parameters
kin_lags     = 3;
kin_rcb_params = struct( ...
    'which_vars', {{'pos','vel','speed'}}, ...
    'rcb_n',      2, ...
    'rcb_hpeaks', [0.05,0.1], ...
    'rcb_b',      0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% array/neuron parameters
in_array     = 'PMd';
out_array    = 'M1';
in_dims      = 16;
out_dims     = 8;
null_size    = ''; %'double',''
pn_kernel_SD = 0.1;
spike_rcb_params = struct( ...
    'which_vars', [in_array '_spikes'], ...
    'rcb_n',      0, ...
    'rcb_hpeaks', [0.01*num_bins,0.1], ...
    'rcb_b',      0.3);
badneuron_params = struct( ...
    'min_fr',1, ...
    'do_shunt_check',1, ...
    'use_trials',{{'epoch','BL'}});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariates parameters
basic_inputs = { ...
    'pos_rcb',3:4; ...
    'vel_rcb',3:4; ...
    'speed_rcb',2};
cov_in.basic = basic_inputs;

% basic_inputs = {};
cov_in.(lower(in_array))  = cat(1,{[in_array '_spikes'],'all'},basic_inputs);
if ~strcmp(in_array,out_array) % potent/null doesn't make sense if they're the same
    cov_in.([lower(in_array) '_pca'])  = cat(1,{[in_array '_pca'],1:in_dims},basic_inputs);
    cov_in.potent  = cat(1,{[in_array out_array '_potent'],'all'},basic_inputs);
    cov_in.null  = cat(1,{[in_array out_array '_null'],'all'},basic_inputs);
end

ANALYSIS_NAME = [in_array out_array '_glm_' name_info];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build function call list for loadTDfiles
func_calls = { ...
    {@getTDidx,'result',use_results}, ...
    {@removeBadTrials,badtrial_params}, ...
    {@removeBadNeurons,badneuron_params}, ...
    {@getSpeed}, ...
    {@binTD,num_bins}, ...
    {@convBasisFunc,kin_rcb_params}, ...
    {@convBasisFunc,spike_rcb_params}};

% get session indices from filedb
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit GLMs
LearningGLM_fit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross validate GLMs
LearningGLM_cv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate GLMs
LearningGLM_eval;


