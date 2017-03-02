% TO DO:
%   1) Predicting in blocks for pR2
%   2) plotting code
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run GLM encoding model
clear;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trial parameters
num_bins  = 5;
badtrial_params = struct( ...
    'ranges', {{'idx_go_cue','idx_movement_on',[5 50]}});
idx_start = {'idx_target_on', 0};
idx_end   = {'idx_trial_end', 0};
train_info = {'BL',[0 0.9]};
test_info = {'BL',[0.9 1]};
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
min_fr    = 5; % in Hz
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
cov_in.full  = cat(1,{[in_array '_spikes'],'all'},cov_in.basic);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
% load('/Users/mattperich/Data/han_td.mat','trial_data');
load('/Users/mattperich/Data/TrialDataFiles/old/Chewie_CO_FF_2016-10-07.mat')
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
td_temp = appendTDs( ...
    truncateAndBin(td,{'idx_target_on',0},{'idx_target_on',50}), ...
    truncateAndBin(td,{'idx_go_cue',0},{'idx_go_cue',50}));
in_dims = estimateDimensionality(td_temp,struct('signal',[in_array '_spikes'],'do_smoothing',true));
out_dims = estimateDimensionality(td_temp,struct('signal',[out_array '_spikes'],'do_smoothing',true));
clear td_temp

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
% Re-bin
td = truncateAndBin(td,num_bins,{'idx_go_cue',0},{'idx_trial_end',0});

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit the models
[glm_info,pr2,rpr2] = deal(struct());
model_names = fieldnames(cov_in);

for m = 1:length(model_names)
    disp(['Calculating ' model_names{m} ' GLM...']); tic;
    [td, glm_info.(model_names{m})] = getModel(td,struct( ...
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
    pr2.(model_names{m}) = evalModel(td,struct( ...
        'model_type','glm', ...
        'out_signals',{glm_info.(model_names{m}).out_signals}, ...
        'model_name',glm_info.(model_names{m}).model_name, ...
        'eval_metric','pr2')); toc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate relative performance of models
temp_names = model_names(~strcmpi(model_names,'basic'));
for m = 1:length(temp_names)
    disp(['Evaluating ' temp_names{m} ' rpR2...']); tic;
    rpr2.(temp_names{m}) = evalModel(td,struct( ...
        'model_type','glm', ...
        'out_signals',{glm_info.(temp_names{m}).out_signals}, ...
        'model_name',{{glm_info.basic.model_name,glm_info.(temp_names{m}).model_name}}, ...
        'eval_metric','pr2')); toc;
end, clear temp_names i m;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some plotting



