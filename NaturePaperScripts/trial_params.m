%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters here so that all paper scripts use the same criteria
% for cutting trials/neurons/etc

%%% for getTDidx('result')
use_results = {'R'};

%%% for removeBadTrialz
badtrial_params = struct();%...
%     'ranges', {{'idx_go_cue','idx_movement_on',[5 50]}});

%%% for removeBadNeurons
badneuron_params = struct( ...
    'min_fr',0.3, ...
    'do_shunt_check',1, ...
    'use_trials',{{'epoch','BL'}});

%%% for trimTD
idx_start = {'idx_target_on', 0};
idx_end   = {'idx_trial_end', -2};

%%% for potent space
pn_kernel_SD = 0.1;

%%% set up initial function calls
trial_func_calls = {...
    {@getTDidx,'result',use_results}, ...
    {@removeBadTrials,badtrial_params}};
