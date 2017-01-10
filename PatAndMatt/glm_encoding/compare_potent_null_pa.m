% load trial_data and params
% plot spiking for a single trial from a single cell as well as predictions
% for that cell from a number of different models
clear;
clc;
close all;

dataSummary;

sessions = { ...
    'Chewie','2016-09-15'; ... % CF
%     'Chewie','2016-10-05'; ...
%     'Chewie','2016-10-07'; ...
%     'Chewie','2016-10-11'; ...
%     'Mihili','2014-02-03'; ...
%     'Mihili','2014-02-17'; ...
%     'Mihili','2014-02-18'; ...
%     'Mihili','2014-03-07'; ...
    };
monkeys = unique(sessions(:,1));
tasks = {'CO'};
perts = {'FF'};
dates = sessions(:,2);

result_codes = {'R'};

neuron_id = 'all';

use_date_idx = ismember(filedb.Monkey,monkeys) & ismember(filedb.Perturbation,perts) & ismember(filedb.Task,tasks);
if ~isempty(dates)
    use_date_idx = use_date_idx & ismember(filedb.Date,dates);
end
use_files = find(use_date_idx);

for idx_file = 1:length(use_files)
    disp(['File ' num2str(idx_file) ' of ' num2str(length(use_files))]);
    % load trial data
    filename = [filedb.Monkey{use_files(idx_file)} '_' filedb.Task{use_files(idx_file)} '_' filedb.Perturbation{use_files(idx_file)} '_' filedb.Date{use_files(idx_file)}];
    load(fullfile(rootDir,TDDir,[filename '.mat']));
    
    load(fullfile(rootDir,TDDir,'trainad_null_move',['FF-PMd-M1_' filename '_cv.mat']),'cv_params');
    params = cv_params;
    params.arrays = {'M1','PMd'};
    [trial_data, params] = glm_process_trial_data(trial_data,params);
    load(fullfile(rootDir,TDDir,'trainad_null_move',['FF-PMd-M1_' filename '.mat']),'params');
    
    
    % get potent/null spaces in baseline
    bl_idx = getTDidx(trial_data,'epoch','bl');
    [V_p1, V_n1, w_cov, w_pred] = getPotentSpace(trial_data(bl_idx),'PMd','M1',params);
    
    % get potent/null spaces in AD
    ad_idx = find(getTDidx(trial_data,'epoch','ad'));
    [V_p2, V_n2, w_cov, w_pred] = getPotentSpace(trial_data(ad_idx(1:floor(0.5*length(ad_idx)))),'PMd','M1',params);
    [V_p3, V_n3, w_cov, w_pred] = getPotentSpace(trial_data(ad_idx(floor(0.5*length(ad_idx))+1:end)),'PMd','M1',params);
    
%     pa = principal_angles(w2(:,1:num_dims),w(:,1:num_dims));
end




