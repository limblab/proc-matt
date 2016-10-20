%% Load data

% Pat's computer
data_fname = 'Mihili_2014-03-03.mat';
data_folder = 'C:\Users\pnlawlor\Box Sync\PatAndMattData';

% Matt's computer

% Load data
load([data_folder '\' data_fname])

%% Load scripts

addpath(genpath('C:\Users\pnlawlor\GoogleDrive\Research\Projects\PatAndMatt\Scripts_PatAndMatt'))

%% Choose visualization parameters

param_struct.trials = 10;


%% Visualize

vis_data(trial_data,param_struct)

