%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run GLM encoding model
clear; clc; close all;
dataSummary;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name for saving
ANALYSIS_NAME = 'PMdM1_glm';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session parameters
monkeys = {'Chewie'};
tasks = {'CO'};
perts = {'FF'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trial parameters
num_bins  = 3;
use_results = {'R','I'};
badtrial_params = struct(...
    'ranges', {{'idx_go_cue','idx_movement_on',[5 50]}});
idx_start = {'idx_target_on', 0};
idx_end   = {'idx_trial_end', 0};
train_info = {'AD',[0.5 1]};
num_boots = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kinematics parameters
kin_rcb_params = struct( ...
    'which_vars', {{'vel'}}, ...
    'rcb_n',      3, ...
    'rcb_hpeaks', [0.01*num_bins,0.2], ...
    'rcb_b',      0.1, ...
    'flip_time',  true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% array/neuron parameters
in_array  = 'PMd';
out_array = 'M1';
badneuron_params = struct( ...
    'min_fr',1, ...
    'do_shunt_check',true);
spike_rcb_params = struct( ...
    'which_vars', [in_array '_spikes'], ...
    'rcb_n',      0, ...
    'rcb_hpeaks', [0.01*num_bins,0.1], ...
    'rcb_b',      0.3);
potentNull_params = struct( ...
    'in_signals',[in_array '_spikes'], ...
    'out_signals',[out_array '_spikes'], ...
    'in_dims',[], ...
    'out_dims',[], ...
    'use_trials',{{'epoch','BL'}}, ...
    'do_smoothing',true, ...
    'trim_idx',{{'idx_movement_on',round(-50/num_bins);'idx_movement_on',round(50/num_bins)}});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM parameters
do_lasso = false;
lasso_lambda = 0.0083;
lasso_alpha = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariates parameters
cov_in.basic = { ...
    'vel','all'; ...
    'vel_rcb','all'};
cov_in.pmd  = cat(1,{[in_array '_pca'],'all'},cov_in.basic);
cov_in.potent  = cat(1,{[in_array out_array '_potent'],'all'},cov_in.basic);
cov_in.null  = cat(1,{[in_array out_array '_null'],'all'},cov_in.basic);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build function call list for loadTDfiles
func_calls = { ...
    {@getTDidx,'result',use_results}, ...
    {@removeBadTrials,badtrial_params}, ...
    {@removeBadNeurons,badneuron_params}, ...
    {@getSpeed}, ...
    {@convBasisFunc,kin_rcb_params}, ...
    {@convBasisFunc,spike_rcb_params}, ...
    {@getPotentSpace,potentNull_params}, ...
    {@trimTD,idx_start,idx_end}, ...
    {@binTD,num_bins}};

session_idx = getFileDBidx(filedb, ...
    {'Task',tasks,'Perturbation',perts,'Monkey',monkeys}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
fnames = cell(1,length(session_idx));
for iFile = 1:length(session_idx)
    file = session_idx(iFile);
    fnames{iFile} = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
end
[trial_data,params] = loadTDfiles(fnames,func_calls{:});


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run GLM
model_names = fieldnames(cov_in);
glm_data = [];
[glm_info,pr2,rpr2] = deal(repmat(struct(),1,length(session_idx)));
for iFile = 1:length(session_idx)
    file = session_idx(iFile);
    
    % get relevant trials
    [~,td] = getTDidx(trial_data, ...
        'date',datestr(filedb.Date{file},'mm-dd-yyyy'), ...
        'monkey',filedb.Monkey{file}, ...
        'task',filedb.Task{file},...
        'perturbation',filedb.Perturbation{file});
    train_trials = getTDidx(td,'epoch',train_info{1},'range',train_info{2});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pull out useful meta data to save with GLM predictions later
    glm_data_file = rmfield(td,setdiff(fieldnames(td),[ ...
        getTDfields(td,'meta'); ...
        getTDfields(td,'unit_guides'); ...
        cellfun(@(x) ['glm_' x],model_names,'uni',0)]));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit the models
    for m = 1:length(model_names)
        disp(['Calculating ' model_names{m} ' GLM...']); tic;
        [td, glm_info(iFile).(model_names{m})] = getModel(td,struct( ...
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
        pr2 = evalModel(td,struct( ...
            'model_type','glm', ...
            'out_signals',{glm_info(iFile).(model_names{m}).out_signals}, ...
            'model_name',glm_info(iFile).(model_names{m}).model_name, ...
            'eval_metric','pr2', ...
            'trial_idx',1:length(td), ...
            'num_boots',num_boots)); toc;
                % add it to the output struct
        for trial = 1:size(pr2,1)
            glm_data_file(trial).([model_names{m} '_pr2']) = squeeze(pr2(trial,:,:));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate relative performance of models
    temp_names = model_names(~strcmpi(model_names,'basic'));
    for m = 1:length(temp_names)
        disp(['Evaluating ' temp_names{m} ' rpR2...']); tic;
        rpr2 = evalModel(td,struct( ...
            'model_type','glm', ...
            'out_signals',{glm_info(iFile).(temp_names{m}).out_signals}, ...
            'model_name',{{glm_info(iFile).basic.model_name,glm_info(iFile).(temp_names{m}).model_name}}, ...
            'eval_metric','pr2', ...
            'trial_idx',1:length(td), ...
            'num_boots',num_boots)); toc;
        % add it to the output struct
        for trial = 1:size(rpr2,1)
            glm_data_file(trial).([temp_names{m} '_rpr2']) = squeeze(rpr2(trial,:,:));
        end
    end, clear temp_names i m;
    
    % put all of the glm predictions together
    glm_data = [glm_data, glm_data_file];
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Package up parameters
params.files = filedb(session_idx,:);
params.glm_info = glm_info;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save results
if ~exist(fullfile(rootDir,resultsDir,ANALYSIS_NAME),'dir')
    mkdir(fullfile(rootDir,resultsDir,ANALYSIS_NAME))
end
save(fullfile(rootDir,resultsDir,ANALYSIS_NAME,[ANALYSIS_NAME '.mat']),'glm_data','params','pr2','rpr2')

%%
clc; close all;

models = {'pmd','null','potent'};

figure; hold all;
for j = 1:length(models)
    a=[];
    for iFile = 1:6
        train_idx = glm_info(iFile).pmd.train_idx;
        
        
        temp = mean(rpr2(iFile).(models{j}),3);
        
        cv_ref = mean(temp(train_idx(1)-30:train_idx(1)-10,:),1);
        idx = cv_ref > 0.01;
        
        temp = moving_average(temp,10);
        
        temp = (temp-repmat(cv_ref,size(temp,1),1))./repmat(cv_ref,size(temp,1),1);
        temp = temp(:,idx);
        a = [a, temp(train_idx(1)-200:train_idx(1)+50,:)];
        
    end
    
    temp = median(a,2);
    plot(temp);
    % plot([train_idx(1) train_idx(1)]-1,[-1 1],'k--');
    % plot([train_idx(end) train_idx(end)]-1,[-1 1],'k--');
    
end


