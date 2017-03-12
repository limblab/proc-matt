%% Prep datasets and get paths etc
clear; clc; close all;
dataSummary;
trial_params;

ANALYSIS_NAME = 'PMdM1_glm_mrt';
file_idx = 1;

load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,'start_params.mat'));

models = {'pmd'};%,'potent','null'};

[~,fname,~] = fileparts(filenames{file_idx});
load(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[fname '_eval.mat']));

train_idx = params.glm_info.(models{m}).train_idx;

