% Use Machens method to esimate dimensionality for M1 and PMd for all data
% files that I have
clear; close all; clc;
dataSummary;

tasks = {'CO'};
perts = {'FF','VR'};

array = 'M1';

session_idx = find( ...
    ismember(filedb.Task,tasks) & ...
    ismember(filedb.Perturbation,perts) & ...
    ~cellfun(@isempty,filedb.FileNames) & ...
    ~ismember(filedb.Monkey,'MrT') & ...
    cellfun(@(x) any(ismember(x,array)),filedb.Arrays));

M1_dims = zeros(1,length(session_idx));
for i = 1:length(session_idx)
    file = session_idx(i);
    disp(['File ' num2str(i) ' of ' num2str(length(session_idx))]);
    
    fname = [filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat'];
    load(fullfile(rootDir,TDDir,fname),'trial_data');
    td = getMoveOnsetAndPeak(trial_data);
    
    [~,td] = getTDidx(td,'result','R','epoch','BL');
    td = removeBadNeurons(td,struct('min_fr',1));
    td = removeBadTrials(td,struct('ranges', ...
        {{'idx_go_cue','idx_movement_on',[0,50]}}));
    
    for j = 1:length(td)
        td(j).target_direction = bin_angles(td(j).target_direction,2*pi/8);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ESTIMATE DIMENSIONALITY
    % td_temp = appendTDs( ...
    %     truncateAndBin(td,{'idx_target_on',0},{'idx_target_on',50}), ...
    %     truncateAndBin(td,{'idx_go_cue',0},{'idx_go_cue',50}));
    td_temp = truncateAndBin(td,{'idx_go_cue',0},{'idx_go_cue',50});
    td_temp = sqrtTransform(td_temp);
    td_temp = smoothSignals(td_temp,struct('signals',[array '_spikes']));
    M1_dims(i) = estimateDimensionality(td_temp,struct('signal',[array '_spikes'],'condition','target_direction'));
    clear td_temp
    
end

%%
array = 'PMd';

session_idx = find( ...
    ismember(filedb.Task,tasks) & ...
    ismember(filedb.Perturbation,perts) & ...
    ~cellfun(@isempty,filedb.FileNames) & ...
    cellfun(@(x) any(ismember(x,array)),filedb.Arrays));

PMd_dims = zeros(1,length(session_idx));
for i = 1:length(session_idx)
    file = session_idx(i);
    disp(['File ' num2str(i) ' of ' num2str(length(session_idx))]);
    
    fname = [filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat'];
    load(fullfile(rootDir,TDDir,fname),'trial_data');
    td = getMoveOnsetAndPeak(trial_data);
    
    [~,td] = getTDidx(td,'result','R','epoch','BL');
    td = removeBadNeurons(td,struct('min_fr',1));
    td = removeBadTrials(td,struct('ranges', ...
        {{'idx_go_cue','idx_movement_on',[0,50]}}));
    
    for j = 1:length(td)
        td(j).target_direction = bin_angles(td(j).target_direction,2*pi/8);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ESTIMATE DIMENSIONALITY
    % td_temp = appendTDs( ...
    %     truncateAndBin(td,{'idx_target_on',0},{'idx_target_on',50}), ...
    %     truncateAndBin(td,{'idx_go_cue',0},{'idx_go_cue',50}));
    td_temp = truncateAndBin(td,{'idx_go_cue',0},{'idx_go_cue',50});
    td_temp = sqrtTransform(td_temp);
    td_temp = smoothSignals(td_temp,struct('signals',[array '_spikes']));
    PMd_dims(i) = estimateDimensionality(td_temp,struct('signal',[array '_spikes'],'condition','target_direction'));
    clear td_temp
    
end

%%
figure;
subplot(121);
hist(M1_dims,1:max([M1_dims,PMd_dims]));
subplot(122);
hist(PMd_dims,1:max([M1_dims,PMd_dims]));