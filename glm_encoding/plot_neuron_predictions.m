% plot spiking for a single trial from a single cell as well as predictions
% for that cell from a number of different models
clear;
clc;
close all;

dataSummary;

monkey = 'Chewie';
date = '2016-10-07';

neuron_id = 'all';

use_file = find(ismember(filedb.Monkey,monkey) & ismember(filedb.Date,date));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the trial data file
filename = [filedb.Monkey{use_file} '_' filedb.Task{use_file} '_' filedb.Perturbation{use_file} '_' filedb.Date{use_file}];
load(fullfile(rootDir,TDDir,[filename '.mat']));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params_file = 'F:\TrialDataFiles\trainad\FF-M1-M1_Chewie_CO_FF_2016-09-15_cv.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the trial data file
load(params_file,'cv_params');
fn = fieldnames(cv_params);
for i = 1:length(fn)
    eval([fn{i} ' = cv_params.' fn{i} ';']);
end, clear i good_cells cov_array pred_array pert epochs pca_w;
params = cv_params;
params.arrays = {'M1','PMd'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[trial_data, params] = glm_process_trial_data(trial_data,params);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop along prediction models and plot trials
cov_array = 'M1';
pred_array = 'M1';

params.pred_array = pred_array;
params.cov_array = cov_array;

if ischar(neuron_id)
    neuron_id = trial_data(1).([pred_array '_unit_guide']);
end



for unit = 35%1:size(neuron_id,1)
    neuron_idx = ismember(trial_data(1).([pred_array '_unit_guide']),neuron_id(unit,:),'rows');
    
    pr2 = zeros(2,2);
    rpr2 = zeros(2,2);
    figure; subplot(1,2,1); hold all;
    load(fullfile(rootDir,TDDir,'trainad',['FF-' cov_array '-' pred_array '_' filename '.mat']),'params');
    [y,x_full,x_basic] = glm_prep_inputs(trial_data,unit,params.train_trials,params.train_start_idx,params.train_end_idx,params);
    b = glmfit(x_full,y,'poisson');
    b_basic = glmfit(x_basic,y,'poisson');
    
    a = cell(size(params.test_epochs,1),1); for i = 1:size(params.test_epochs,1), a{i} = params.test_epochs(i,:); end
    params.test_epochs = a; clear a i;
    trial_idx = find(strcmpi(params.test_epochs,'AD'),30,'first')';
    trial_idx = trial_idx(5:end);
    
    [y,x_full,x_basic] = glm_prep_inputs(trial_data,unit,trial_idx,params.train_start_idx,params.train_end_idx,params);
    plot(y,'k-','LineWidth',2);
    yfit = exp([ones(size(x_full,1),1), x_full]*b)';
    yfit_basic = exp([ones(size(x_basic,1),1), x_basic]*b_basic)';
    plot(yfit_basic,'LineWidth',2);
    plot(yfit,'LineWidth',2);
    bs = randi(length(yfit),length(yfit),num_bootstraps);
    pr2(1,:) = prctile(compute_pseudo_R2(y(bs),yfit(bs),mean(y)),[2.5 97.5]);
    rpr2(1,:) = prctile(compute_rel_pseudo_R2(y(bs),yfit_basic(bs),yfit(bs)),[2.5 97.5]);
    
    legend({'Real',['Basic, pr2 = ' num2str(mean(pr2(1,:)),3)],['M1-M1, rpr2 = ' num2str(mean(rpr2(1,:)),3)]});
    set(gca,'Box','off','TickDir','out','FontSize',14); axis('tight');
    ylabel('Spike Count','FontSize',14);
    xlabel('Time bins','FontSize',14);
    title('First 10 CF Trials','FontSize',14);
    
    if mean(rpr2(1,:)) > 0.3
        
        subplot(1,2,2); hold all;
        load(fullfile(rootDir,TDDir,'trainad',['FF-' cov_array '-' pred_array '_' filename '.mat']),'params');
        [y,x_full,x_basic] = glm_prep_inputs(trial_data,unit,params.train_trials,params.train_start_idx,params.train_end_idx,params);
        b = glmfit(x_full,y,'poisson');
        b_basic = glmfit(x_basic,y,'poisson');
        
        a = cell(size(params.test_epochs,1),1); for i = 1:size(params.test_epochs,1), a{i} = params.test_epochs(i,:); end
        params.test_epochs = a; clear a i;
        trial_idx = find(strcmpi(params.test_epochs,'AD'),80,'first')';
        trial_idx = trial_idx(50:end);
        
        [y,x_full,x_basic] = glm_prep_inputs(trial_data,unit,trial_idx,params.train_start_idx,params.train_end_idx,params);
        plot(y,'k-','LineWidth',2);
        yfit = exp([ones(size(x_full,1),1), x_full]*b)';
        yfit_basic = exp([ones(size(x_basic,1),1), x_basic]*b_basic)';
        plot(yfit_basic,'LineWidth',2);
        plot(yfit,'LineWidth',2);
        bs = randi(length(yfit),length(yfit),num_bootstraps);
        pr2(2,:) = prctile(compute_pseudo_R2(y(bs),yfit(bs),mean(y)),[2.5 97.5]);        
        rpr2(2,:) = prctile(compute_rel_pseudo_R2(y(bs),yfit_basic(bs),yfit(bs)),[2.5 97.5]);
        
        legend({'Real',['Basic, pr2 = ' num2str(mean(pr2(2,:)),3)],['M1-M1, rpr2 = ' num2str(mean(rpr2(2,:)),3)]});
        set(gca,'Box','off','TickDir','out','FontSize',14); axis('tight');
        ylabel('Spike Count','FontSize',14);
        xlabel('Time bins','FontSize',14);
        title('Last 10 CF Trials','FontSize',14);
        
        pause;
    end
    close all;
    
end



%% Compare potent/null
for unit = 1:size(neuron_id,1)
    neuron_idx = ismember(trial_data(1).([pred_array '_unit_guide']),neuron_id(unit,:),'rows');
    
    pr2 = zeros(2,2);
    figure; subplot(1,2,1); hold all;
    load(fullfile(rootDir,TDDir,'trainad_potent',['FF-PMd-M1_' filename '.mat']),'params');
    [y,x_full,x_basic] = glm_prep_inputs(trial_data,unit,params.train_trials,params.train_start_idx,params.train_end_idx,params);
    b = glmfit(x_full,y,'poisson');
    
    a = cell(size(params.test_epochs,1),1); for i = 1:size(params.test_epochs,1), a{i} = params.test_epochs(i,:); end
    params.test_epochs = a; clear a i;
    trial_idx = find(strcmpi(params.test_epochs,'AD'),10,'first')';
    
    [y,x_full,~] = glm_prep_inputs(trial_data,unit,trial_idx,params.train_start_idx,params.train_end_idx,params);
    plot(y,'k-','LineWidth',2);
    yfit = exp([ones(size(x_full,1),1), x_full]*b)';
    plot(yfit,'LineWidth',2);
    bs = randi(length(yfit),length(yfit),num_bootstraps);
    pr2(1,:) = prctile(compute_pseudo_R2(y(bs),yfit(bs),mean(y)),[2.5 97.5]);
    
    if mean(pr2(1,:)) > 0.1
        load(fullfile(rootDir,TDDir,'trainad_null',['FF-PMd-M1_' filename '.mat']),'params');
        [y,x_full,~] = glm_prep_inputs(trial_data,unit,params.train_trials,params.train_start_idx,params.train_end_idx,params);
        b = glmfit(x_full,y,'poisson');
        [y,x_full,~] = glm_prep_inputs(trial_data,unit,trial_idx,params.train_start_idx,params.train_end_idx,params);
        yfit = exp([ones(size(x_full,1),1), x_full]*b)';
        plot(yfit,'LineWidth',2);
        bs = randi(length(yfit),length(yfit),num_bootstraps);
        pr2(2,:) = prctile(compute_pseudo_R2(y(bs),yfit(bs),mean(y)),[2.5 97.5]);
        legend({'Real',['potent: ' num2str(mean(pr2(1,:)),3)],['null: ' num2str(mean(pr2(2,:)),3)]});
        set(gca,'Box','off','TickDir','out','FontSize',14); axis('tight');
        ylabel('Spike Count','FontSize',14);
        xlabel('Time bins','FontSize',14);
        title('First 10 CF Trials','FontSize',14);
        
        subplot(1,2,2); hold all;
        load(fullfile(rootDir,TDDir,'trainad_potent',['FF-PMd-M1_' filename '.mat']),'params');
        [y,x_full,~] = glm_prep_inputs(trial_data,unit,params.train_trials,params.train_start_idx,params.train_end_idx,params);
        b = glmfit(x_full,y,'poisson');
        
        a = cell(size(params.test_epochs,1),1); for i = 1:size(params.test_epochs,1), a{i} = params.test_epochs(i,:); end
        params.test_epochs = a; clear a i;
        trial_idx = find(strcmpi(params.test_epochs,'AD'),10,'last')';
        
        [y,x_full,~] = glm_prep_inputs(trial_data,unit,trial_idx,params.train_start_idx,params.train_end_idx,params);
        plot(y,'k-','LineWidth',2);
        yfit = exp([ones(size(x_full,1),1), x_full]*b)';
        plot(yfit,'LineWidth',2);
        bs = randi(length(yfit),length(yfit),num_bootstraps);
        pr2(1,:) = prctile(compute_pseudo_R2(y(bs),yfit(bs),mean(y)),[2.5 97.5]);        
        
        load(fullfile(rootDir,TDDir,'trainad_null',['FF-PMd-M1_' filename '.mat']),'params');
        [y,x_full,x_basic] = glm_prep_inputs(trial_data,unit,params.train_trials,params.train_start_idx,params.train_end_idx,params);
        b = glmfit(x_full,y,'poisson');
        [y,x_full,~] = glm_prep_inputs(trial_data,unit,trial_idx,params.train_start_idx,params.train_end_idx,params);
        yfit = exp([ones(size(x_full,1),1), x_full]*b)';
        plot(yfit,'LineWidth',2);
        bs = randi(length(yfit),length(yfit),num_bootstraps);
        pr2(2,:) = prctile(compute_pseudo_R2(y(bs),yfit(bs),mean(y)),[2.5 97.5]);
        
        legend({'Real',['potent: ' num2str(mean(pr2(1,:)),3)],['null: ' num2str(mean(pr2(2,:)),3)]});
        set(gca,'Box','off','TickDir','out','FontSize',14); axis('tight');
        ylabel('Spike Count','FontSize',14);
        xlabel('Time bins','FontSize',14);
        title('Last 10 CF Trials','FontSize',14);
        
        pause;
    end
    close all;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




