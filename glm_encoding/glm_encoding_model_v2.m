%%
% STUFF TO DO
% x - PCA as input
% x - null and potent space as input
% x - mahalonobis distance at time (e.g. target presentation, go cue, peak
%   speed, etc) between baseline and force as function of adaptation (for
%   PCA of M1, PMd, etc)
%   - compare cell tuning. Do predictions with only "untuned" cells as
%   covariates and show that they have information
%   - include M1 and PMd, look at what PMd adds over M1
% x - include target direction as covariate with kinematics
%   - are difficult to predict neurons "soloist" cells? Do R2 as function of norm2 of the PC weights
%   - normalize by variance in some way instead of CV mean
%   - correlate predictions of null/potent space instead of R2?
%   - Do predictions as function of M1/PMd dimensions
%   - Dig up old Mihili file without washout but with good CF
%   - add shunt check to trial_data processing code
%   - predict potent space and null space trial by trial
%   - predict just movement
%
% FUTURE ANALYSES
%   - Correlate errors with projections into null space?
%   - look at relationship between null/potent weights and tuning?
%   - Look for patterns in correlograms while sorting. See if this predicts
%   any behavior in interactions between neurons
%   - look at spatial location on array with manifold


%%% Differences of current trial_data from old:
%       Old one cut trials < 0.1 ms
%       Old one cut trials where the movement or peak weren't found
%

clear;
clc;
close all;

dataSummary;

outputSubdir = 'trainad_potent';
params_comment = 'training on AD data with 50 msec bins';

params_file = '';% optionally give .mat file containing params to use
% Must be in the outputSubdir

%       'Chewie','2016-09-19'; ...
          
sessions = { ...
    'Chewie','2016-09-15'; ... % CF
    'Chewie','2016-09-21'; ...
    'Chewie','2016-10-05'; ...
    'Chewie','2016-10-07'; ...
    'Chewie','2016-10-11'; ...
    'Mihili','2014-02-03'; ...
    'Mihili','2014-02-17'; ...
    'Mihili','2014-02-18'; ...
    'Mihili','2014-03-07'; ...
    'Chewie','2016-09-09'; ... % VR
    'Chewie','2016-09-12'; ...
    'Chewie','2016-09-14'; ...
    'Chewie','2016-10-06'; ...
    'Mihili','2014-03-03'; ...
    'Mihili','2014-03-04'; ...
    'Mihili','2014-03-06'; ...
    };


monkeys = unique(sessions(:,1));
tasks = {'CO'};
perts = {'FF'};
dates = sessions(:,2);
% {COVARIATE, PREDICTED}, will loop along rows
% array_pairs = {'PMd','M1';'M1','M1';'PMd','PMd'};
array_pairs = {'PMd','M1'};

%%
if isempty(params_file)
    dt                     = 0.01; % time step size for data
    bin_size               = 5;    % how many samples to group together when rebinning (bin width is num_samples*dt)
    
    result_codes           = {'R'}; % which trials to include
    
    training_data          = {'AD',[0.5 1]}; % either proportion range e.g. [0 0.5] or # trails from start (100) or end (-100)
    %     testing_data           = {'BL','all'; ... % should match this to be proportion or trial range based on training_data
    %         'AD',[0 0.5]; ...
    %         'WO','all'};
    testing_data           = {'AD',[0 0.5]};
    block_size_testing     = 1; % size of blocks in # trials for AD/WO
    % note: can change this and still reload CV
    
    % how to truncate trials {idx name, number of bins after}
    train_start_idx        = {'idx_target_on',0}; %{'idx_target_on',-4}
    train_end_idx          = {'idx_trial_end',-2}; %{'idx_go_cue',4}
    test_start_idx         = {'idx_target_on',0}; % NOTE: CAN CHANGE THESE
    test_end_idx           = {'idx_trial_end',-2}; %   AND STILL RELOAD CV
    %   NOTE: this is after rebinning at the moment
    
    cv_folds               = 10; % how many folds for cross-validation
    cv_rand                = true; % whether to randomly select CV datapoints
    
    block_size_fr_test     = 100;   % how many trials to group for FR test
    fr_test_alpha          = 1e-3; % p value cut off for t-test
    fr_min                 = 0.15; % minimum session-wide spiking for inclusion
    
    do_lasso               = false;  % if you want to regularize
    lasso_lambda           = 0.0083; % lambda parameter for lasso
    lasso_alpha            = 0.01;   % alpha parameter for lasso
    
    do_pca                 = 'potent'; % false, 'pca', 'potent', 'null'
    smooth_pca_proj_spikes = false;
    center_pca             = true;
    pca_dims               = struct('M1',1:8,'PMd',1:16); % ...'ARRAY','all' or specify which dimensions, e.g. 1:30...
    potent_trials          = {'epoch','bl'};
    
    do_rcb                 = true;  % use raised cosine basis functions for history
    do_all_history         = false; % include history for all units
    do_self_history        = false;
    unit_lags              = 2;     % how many bins in the past for each neuron to include as covariates (if ~do_rcb)
    rcb_hpeaks             = [dt*bin_size,0.1];
    rcb_b                  = 0.1;
    %   NOTE: this is after rebinning at the moment
    
    kin_signals            = {'pos','vel','speed'}; % names of kinematic covariates
    kin_lags               = 3; % how many bins to lag (if ~do_rcb) or how many raised cosines (if do_rcb)
    filter_velocity        = false;
    %   NOTE: this is after rebinning at the moment
    
    num_bootstraps         = 1000; % how many resamples for R2
    
else % load up params from whatever file
    load(params_file,'cv_params');
    fn = fieldnames(params);
    for i = 1:length(fn)
        eval([fn{i} ' = cvparams.' fn{i} ';']);
    end, clear i good_cells cov_array pred_array pert epochs pca_w;
end



%% Here is where I check all of the inputs
possible_arrays = {'M1','PMd','M1PMd'}; % list all here!

if bin_size < 1 || block_size_testing < 1
    error('Invalid bin size or testing block size');
elseif bin_size == 1
    disp('|-----------------------------------------------------------|')
    disp('| BIN SIZE IS 1... SOMETIMES THIS CRASHES THINGS. CONTINUE? |');
    disp('|-----------------------------------------------------------|')
    pause;
end

if ischar(do_pca) && ~all(ismember(fieldnames(pca_dims),possible_arrays))
    error('PCA Dimension array name not recognized...');
end

if unit_lags < 0, error('What are we supposed to do with a negative lag?'); end

if ~exist('params_comment','var'), params_comment = ''; end


%%
if ~exist(fullfile(rootDir,resultsDir,outputSubdir),'dir')
    mkdir(fullfile(rootDir,resultsDir,outputSubdir));
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% disp(['Training on ' training_epoch{1} '...']);

for idx_pert = 1:length(perts)
    pert = perts{idx_pert};
    disp(['Starting perturbation ' pert '...']);
    
    use_date_idx = ismember(filedb.Monkey,monkeys) & ismember(filedb.Perturbation,perts(idx_pert)) & ismember(filedb.Task,tasks);
    if ~isempty(dates)
        use_date_idx = use_date_idx & ismember(filedb.Date,dates);
    end
    use_files = find(use_date_idx);
    
    for idx_arrays = 1:size(array_pairs,1)
        
        disp(['STARTING ARRAY PAIR ' array_pairs{idx_arrays,1} ' -> ' array_pairs{idx_arrays,2} '...']);
        
        cov_array = array_pairs{idx_arrays,1};
        pred_array = array_pairs{idx_arrays,2};
        
        for idx_file = 1:length(use_files)
            disp(['File ' num2str(idx_file) ' of ' num2str(length(use_files)) '...'])
            %             close all;
            %             figure;
            
            epochs = filedb.Epochs{use_files(idx_file)};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Package up parameters for saving
            clear params;
            params.cov_array = cov_array;
            params.pred_array = pred_array;
            params.training_data = training_data;
            params.testing_data = testing_data;
            params.block_size_fr_test = block_size_fr_test;
            params.fr_test_alpha = fr_test_alpha;
            params.fr_min = fr_min;
            params.dt = dt;
            params.bin_size = bin_size;
            params.cv_folds = cv_folds;
            params.cv_rand = cv_rand;
            params.do_rcb = do_rcb;
            params.do_self_history = do_self_history;
            params.do_all_history = do_all_history;
            if do_rcb
                params.rcb_hpeaks = rcb_hpeaks;
                params.rcb_b = rcb_b;
            end
            params.unit_lags = unit_lags;
            params.kin_lags = kin_lags;
            params.do_pca = do_pca;
            params.smooth_pca_proj_spikes = smooth_pca_proj_spikes;
            params.center_pca = center_pca;
            params.potent_trials = potent_trials;
            if ischar(do_pca)
                params.pca_dims = pca_dims;
            end
            params.num_bootstraps = num_bootstraps;
            params.do_lasso = do_lasso;
            if do_lasso
                params.lasso_lambda = lasso_lambda;
                params.lasso_alpha = lasso_alpha;
            end
            params.result_codes = result_codes;
            params.kin_signals = kin_signals;
            params.filter_velocity = filter_velocity;
            params.train_start_idx = train_start_idx;
            params.train_end_idx = train_end_idx;
            cv_params = params; % these will be saved separately
            params.arrays = unique(array_pairs(idx_arrays,:));
            params.pert = pert;
            params.filedb = filedb;
            params.epochs = epochs;
            params.test_start_idx = test_start_idx;
            params.test_end_idx = test_end_idx;
            params.block_size_testing = block_size_testing;
            params.comment = params_comment;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if params.smooth_pca_proj_spikes
                disp('SMOOTHING PCA PROJECTION SPIKES!');
            end
            if params.center_pca
                disp('CENTERING PCA!');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load the trial data file
            filename = [filedb.Monkey{use_files(idx_file)} '_' filedb.Task{use_files(idx_file)} '_' filedb.Perturbation{use_files(idx_file)} '_' filedb.Date{use_files(idx_file)}];
            load(fullfile(rootDir,TDDir,[filename '.mat']));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Filter velocity traces at 14 Hz
            if filter_velocity
            if strcmpi(pert,'FF'), disp('FILTERING VELOCITY!'); end
            for trial = 1:length(trial_data)
                if strcmpi(trial_data(trial).epoch,'AD') && strcmpi(pert,'FF')
                    fs = 100;
                    Wn = 14*2/fs;
                    deg = 3;
                    [B,A] = butter(deg,Wn,'low');
                    trial_data(trial).vel(:,1) = filter(B,A,trial_data(trial).vel(:,1));
                    trial_data(trial).vel(:,2) = filter(B,A,trial_data(trial).vel(:,2));
                    % integrate to get new position
%                     trial_data(trial).pos(:,1) = cumtrapz(trial_data(trial).vel(:,1));
%                     trial_data(trial).pos(:,2) = cumtrapz(trial_data(trial).vel(:,2));
                end
            end
            end
            
            trimming_weird_trials;
            

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % will save cross validated BL predictions so I don't have to
            % redo them to to apply the model to other time points
            cv_fn = fullfile(rootDir,resultsDir,outputSubdir,[pert '-' cov_array '-' pred_array '_' filename '_cv.mat']);
            do_cv = true;
            if exist(cv_fn,'file')
                temp_params = cv_params;
                load(cv_fn,'cv_params');
                % check all parameters. Must be exactly the same.
                [~,d,~]=comp_struct(temp_params,cv_params);
                cv_params = temp_params;
                if isempty(d) || isempty(fieldnames(d))
                    disp('Compatible CV file found. Reloading.');
                    do_cv = false;
                    load(cv_fn,'pr2_full_cv','pr2_basic_cv','rpr2_cv');
                end, clear d params2;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmpi(filedb.Monkey{use_files(idx_file)},'Mihili')
                [~,trial_data] = getTDidx(trial_data,'result',result_codes);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % bin and prepare signals
            [trial_data, params] = glm_process_trial_data(trial_data,params);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % filter trials and partition testing/training sets
            [~,trial_data] = getTDidx(trial_data,'result',result_codes);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Duplicating trial data for test set')
            trial_data_test = trial_data;
            trial_data_test = rmfield(trial_data_test,getTDfields(trial_data_test,'time'));
            
            
            if length(training_data{2}) > 1 % use percentages
                [train_trials,test_trials, test_epochs] = deal([]);
                for iBlock = 1:size(training_data,1)
                    temp = find(strcmpi({trial_data.epoch},training_data{iBlock,1}));
                    if ischar(training_data{iBlock,2}) % use all
                        train_trials = [train_trials, temp];
                    else
                        train_trials = [train_trials, temp( 1+floor(length(temp)*training_data{iBlock,2}(1)):floor(length(temp)*training_data{iBlock,2}(2)) )];
                    end
                end
                for iBlock = 1:size(testing_data,1)
                    temp = find(strcmpi({trial_data.epoch},testing_data{iBlock,1}));
                    if ~ischar(testing_data{iBlock,2}) % use all
                        temp = temp( 1+floor(length(temp)*testing_data{iBlock,2}(1)):floor(length(temp)*testing_data{iBlock,2}(2)) );
                    end
                    temp = reshape(temp(1:end-rem(length(temp),block_size_testing)),block_size_testing,(length(temp) - rem(length(temp),block_size_testing))/block_size_testing)';
                    test_trials = [test_trials; temp];
                    test_epochs = [test_epochs; repmat(testing_data{iBlock,1},size(temp,1),1)];
                end
            else % use number of trials
                [train_trials,test_trials, test_epochs] = deal([]);
                for iBlock = 1:size(testing_data,1)
                    temp = find(strcmpi({trial_data.epoch},testing_data{iBlock,1}));
                    if ~ischar(testing_data{iBlock,2}) % use all
                        temp = temp( testing_data{iBlock,2}(1):testing_data{iBlock,2}(2) );
                    end
                    temp = reshape(temp(1:end-rem(length(temp),block_size_testing)),block_size_testing,(length(temp) - rem(length(temp),block_size_testing))/block_size_testing)';
                    test_trials = [test_trials; temp];
                    test_epochs = [test_epochs; repmat(testing_data{iBlock,1},size(temp,1),1)];
                end
                for iBlock = 1:size(training_data,1)
                    temp = find(strcmpi({trial_data.epoch},training_data{iBlock,1}));
                    if ischar(training_data{iBlock,2}) % use all
                        train_trials = [train_trials, temp];
                    else
                        temp_idx = false(1,length(temp));
                        if any(strcmpi(testing_data(:,1),training_data{1,1}))
                            test_idx = testing_data{strcmpi(testing_data(:,1),training_data{1,1}),2}(1):testing_data{strcmpi(testing_data(:,1),training_data{1,1}),2}(2);
                        else
                            test_idx = false(1,length(temp));
                        end
                        
                        if training_data{2} < 0
                            train_idx = length(temp)+training_data{2}+1:length(temp);
                        else
                            train_idx = 1:training_data{2};
                        end
                        
                        temp_idx(train_idx) = true;
                        
                        if any(ismember(test_idx, train_idx));
                            disp('Some overlap in training and testing...');
                            temp_idx(test_idx) = false;
                        end
                        
                        % in the event there are not enough trials for the testing and training sets, make sure they don't overlap
                        train_trials = [train_trials, temp( temp_idx )];
                    end
                end
            end
            params.train_trials = train_trials;
            params.test_trials = test_trials;
            params.test_epochs = test_epochs;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % At last the data is processed!
            num_pred_neurons = size(trial_data(1).([pred_array '_spikes']),2);
            num_cov_neurons = size(trial_data(1).([cov_array '_spikes']),2);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute PCA if we don't need to leave out the pred neuron
            if ~strcmpi(cov_array,pred_array) && ischar(do_pca) % we don't have to leave one out so just fit once
                if strcmpi(do_pca,'pca')
                    [~,temp] = getPCA(trial_data,struct('signals',{{[cov_array '_spikes'],1:num_cov_neurons}},'use_trials',1:length(trial_data)));
                    params.pca_w = temp.w;
                    params.pca_mu = temp.mu;
                elseif strcmpi(do_pca,'potent') || strcmpi(do_pca,'null')
                    disp('Getting potent/null space.')
                    % check that PMd is cov_array
                    if ~strcmpi(cov_array,'PMd') || ~strcmpi(pred_array,'M1')
                        error('PMd must be covariate array and M1 must be predicted array for potent/null space analysis');
                    end
                    % do SVD to get null/potent spaces
                    pn_params = struct( ...
                        'in_signals', [cov_array '_spikes'], ...
                        'out_signals',[pred_array '_spikes'], ...
                        'in_dims',max(params.pca_dims.(cov_array)), ...
                        'out_dims',max(params.pca_dims.(pred_array)), ...
                        'do_smoothing',false, ...
                        'use_trials',{potent_trials});
                    
                    [~,temp] = getPotentSpace(trial_data,pn_params);
                    params.pca_w = temp.w_in;
                    params.pca_mu = temp.mu_in;
                    params.V_potent = temp.V_potent;
                    params.V_null = temp.V_null;
                end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Loop along neurons
            if do_cv
                disp('No CV file found or parameters did not match. Running CV...');
                [pr2_full_cv,rpr2_cv,pr2_basic_cv] = deal(NaN(num_pred_neurons,cv_folds,2));
            end
            bl_unit_models = repmat(struct('b_basic',[],'d_basic',[],'s_basic',[],'b_full',[],'d_full',[],'s_full',[],'ymean',[]),num_pred_neurons,1);
            [pr2_full,rpr2,pr2_basic] = deal(NaN(num_pred_neurons,size(test_trials,1),2));
            for unit = 1:num_pred_neurons
                tic;
                disp(['Running neuron ' num2str(unit) ' of ' num2str(num_pred_neurons)]);
                
                if strcmpi(cov_array,pred_array) && strcmpi(do_pca,'pca')  % recompute PCA each time
                    % leave predicted neuron out of PCA
                    neuron_idx = 1:num_pred_neurons; neuron_idx(unit) = [];
                    [~,temp] = getPCA(trial_data,struct('signals',{{[cov_array '_spikes'],neuron_idx}},'use_trials',1:length(trial_data)));
                    % used in glm_prep_inputs
                    params.pca_w = temp.w;
                    params.pca_mu = temp.mu;
                    % note that if it's the same array, null/potent doesn't make sense
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Build baseline model
                [y,x_full,x_basic] = glm_prep_inputs(trial_data,unit,train_trials,train_start_idx,train_end_idx,params);
                
                if do_lasso
                    [b_basic,d_basic,s_basic,d_full,s_full] = deal([]);
                    [b,s] = lassoglm(x_full,y,'poisson','lambda',lasso_lambda,'alpha',lasso_alpha);
                    b_full = [s.Intercept; b];
                    if ~isempty(x_basic)
                        [b,s] = lassoglm(x_basic,y,'poisson','lambda',lasso_lambda,'alpha',lasso_alpha);
                        b_basic = [s.Intercept; b];
                    end
                else
                    
                    [b_full,d_full,s_full] = glmfit(x_full,y,'poisson');
                    
                    if ~isempty(x_basic)
                        [b_basic,d_basic,s_basic] = glmfit(x_basic,y,'poisson');
                    else
                        [b_basic,d_basic,s_basic] = deal([]);
                    end
                end
                
                bl_unit_models(unit).b_basic = b_basic;
                bl_unit_models(unit).d_basic = d_basic;
                bl_unit_models(unit).s_basic = s_basic;
                bl_unit_models(unit).b_full = b_full;
                bl_unit_models(unit).d_full = d_full;
                bl_unit_models(unit).s_full = s_full;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Cross validate
                if do_cv
                    
                    if cv_rand
                        folds = reshape(randperm(size(y,1)-rem(size(y,1),cv_folds)),[],cv_folds)';
                    else
                        folds = reshape(1:size(y,1)-rem(size(y,1),cv_folds),[],cv_folds)';
                    end
                    
                    for fold = 1:cv_folds
                        test_idx = folds(fold,:);
                        train_idx = true(1,size(y,1)); train_idx(test_idx) = false;
                        
                        y_train = y(train_idx);
                        if do_lasso
                            [b,s] = lassoglm(x_full(train_idx,:),y_train,'poisson','lambda',lasso_lambda,'alpha',lasso_alpha);
                            b_full = [s.Intercept; b];
                            if ~isempty(b_basic)
                                [b,s] = lassoglm(x_basic(train_idx,:),y_train,'poisson','lambda',lasso_lambda,'alpha',lasso_alpha);
                                b_basic = [s.Intercept; b];
                            end
                        else
                            b_full = glmfit(x_full(train_idx,:),y_train,'poisson');
                            if ~isempty(b_basic)
                                b_basic = glmfit(x_basic(train_idx,:),y_train,'poisson');
                            end
                        end
                        
                        y_test = y(test_idx);
                        yfit_full = exp([ones(length(test_idx),1), x_full(test_idx,:)]*b_full)';
                        if ~isempty(b_basic)
                            yfit_basic = exp([ones(length(test_idx),1), x_basic(test_idx,:)]*b_basic)';
                        end
                        
                        % compute R2
                        % get vector of bootstrap indices
                        bs = randi(length(y_test),length(y_test),num_bootstraps);
                        pr2_full_cv(unit,fold,:) = prctile(compute_pseudo_R2(y_test(bs),yfit_full(bs),mean(y_test)),[2.5 97.5]);
                        if ~isempty(b_basic)
                            rpr2_cv(unit,fold,:) = prctile(compute_rel_pseudo_R2(y_test(bs),yfit_basic(bs),yfit_full(bs)),[2.5 97.5]);
                            pr2_basic_cv(unit,fold,:) = prctile(compute_pseudo_R2(y_test(bs),yfit_basic(bs),mean(y_test)),[2.5 97.5]);
                        end
                    end
                    disp(['pR2 = ' num2str(mean(mean(pr2_full_cv(unit,:,:),3)))]);
                    disp(['rpR2 = ' num2str(mean(mean(rpr2_cv(unit,:,:),3)))]);
                    
                    if 0
                        figure;
                        hold all;
                        plot(y_test,'LineWidth',2);
                        plot(yfit_full,'LineWidth',2);
                        %plot(yfit_basic);
                        title(['pseudo-r2 = ' num2str(mean(mean(pr2_full_cv(unit,:,:),3)))],'FontSize',14);
                        set(gca,'Box','off','TickDir','out','FontSize',14);
                        %legend({['rpr2=' num2str(mean(rpr2_cv(unit,end,:)))],['pr2_f=' num2str(mean(pr2_full_cv(unit,end,:)))],['pr2_b=' num2str(mean(pr2_basic_cv(unit,end,:)))]});
                        pause;
                        close all;
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Test throughout session on new data
                for e = 1:size(test_trials,1)
                    %disp(['Predicting ' num2str(e) ' of ' num2str(size(testing_idx,1)) '...']);
                    [y,x_full,x_basic] = glm_prep_inputs(trial_data,unit,test_trials(e,:),test_start_idx,test_end_idx,params);
                    
                    yfit_full = exp([ones(size(x_full,1),1), x_full]*bl_unit_models(unit).b_full)';
                    if ~isempty(b_basic)
                        yfit_basic = exp([ones(size(x_basic,1),1), x_basic]*bl_unit_models(unit).b_basic)';
                    end
                    
                    if block_size_testing == 1
                        % add prediction to trial_data
                        trial_data_test(test_trials(e,:)).y(:,unit) = y;
                        trial_data_test(test_trials(e,:)).x_full = x_full;
                        trial_data_test(test_trials(e,:)).yfit_full(:,unit) = yfit_full;
                        if ~isempty(b_basic)
                            trial_data_test(test_trials(e,:)).x_basic = x_basic;
                            trial_data_test(test_trials(e,:)).yfit_basic(:,unit) = yfit_basic;
                        end
                    end
                    
                    % compute R2
                    % get vector of bootstrap indices
                    bs = randi(length(y),length(y),num_bootstraps);
                    pr2_full(unit,e,:) = prctile(compute_pseudo_R2(y(bs),yfit_full(bs),mean(y)),[2.5 97.5]);
                    if ~isempty(b_basic)
                        rpr2(unit,e,:) = prctile(compute_rel_pseudo_R2(y(bs),yfit_basic(bs),yfit_full(bs)),[2.5 97.5]);
                        pr2_basic(unit,e,:) = prctile(compute_pseudo_R2(y(bs),yfit_basic(bs),mean(y)),[2.5 97.5]);
                    end

                end
                toc;
            end % end neuron loop
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Save results
            if do_cv
                save(fullfile(rootDir,resultsDir,outputSubdir,[pert '-' cov_array '-' pred_array '_' filename '_cv.mat']),'pr2_full_cv','pr2_basic_cv','rpr2_cv','cv_params');
            end
            
            results = struct( ...
                'pr2_full', pr2_full, ...
                'pr2_basic', pr2_basic, ...
                'rpr2', rpr2, ...
                'pr2_full_cv', pr2_full_cv, ...
                'pr2_basic_cv', pr2_basic_cv, ...
                'rpr2_cv', rpr2_cv, ...
                'bl_model', bl_unit_models);
            
            save(fullfile(rootDir,resultsDir,outputSubdir,[pert '-' cov_array '-' pred_array '_' filename '_td.mat']),'trial_data_test','params');
            save(fullfile(rootDir,resultsDir,outputSubdir,[pert '-' cov_array '-' pred_array '_' filename '.mat']),'results','params');
            clear file_results;
        end % end array pair loop
    end % end file loop
end % end perturbation loop
