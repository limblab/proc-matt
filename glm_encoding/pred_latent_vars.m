%%
% Predict potent space from null space


clear;
clc;
close all;

dataSummary;

outputSubdir = 'trainad_potent_null_pred';
params_comment = 'predicting potent space from null space during learning';

params_file = '';% optionally give .mat file containing params to use
% Must be in the outputSubdir

sessions = { ...
    %     'Chewie','2016-09-09'; ... % VR
    %     'Chewie','2016-09-12'; ...
    %     'Chewie','2016-09-14'; ...
    %     'Chewie','2016-10-06'; ...
    %     'Mihili','2014-03-03'; ...
    %     'Mihili','2014-03-04'; ...
    %     'Mihili','2014-03-06'; ...
%         'Chewie','2016-09-15'; ... % CF
%         'Chewie','2016-10-05'; ...
%     'Chewie','2016-10-07'; ...
%         'Chewie','2016-10-11'; ...
        'Mihili','2014-02-03'; ...
        'Mihili','2014-02-17'; ...
        'Mihili','2014-02-18'; ...
        'Mihili','2014-03-07'; ...
    };


monkeys = unique(sessions(:,1));
tasks = {'CO'};
perts = {'FF','VR'};
dates = sessions(:,2);
% {COVARIATE, PREDICTED}, will loop along rows
cov_array = 'PMd';
pred_array = 'M1';

%%
if isempty(params_file)
    dt                     = 0.01; % time step size for data
    bin_size               = 5;    % how many samples to group together when rebinning (bin width is num_samples*dt)
    
    result_codes           = {'R'}; % which trials to include
    
    training_data          = {'AD',[0.5 1]}; % either proportion range e.g. [0 0.5] or # trails from start (100) or end (-100)
    testing_data           = {'BL','all'
        'AD',[0 0.5]; ... % should match this to be proportion or trial range based on training_data
        'WO','all'};
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
    
    smooth_pca_proj_spikes = true;
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
    
    kin_signals            = {'pos','vel','speed','force'}; % names of kinematic covariates
    kin_lags               = 3; % how many bins to lag (if ~do_rcb) or how many raised cosines (if do_rcb)
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
        params.smooth_pca_proj_spikes = smooth_pca_proj_spikes;
        params.center_pca = center_pca;
        params.potent_trials = potent_trials;
        params.pca_dims = pca_dims;
        params.num_bootstraps = num_bootstraps;
        params.do_lasso = do_lasso;
        if do_lasso
            params.lasso_lambda = lasso_lambda;
            params.lasso_alpha = lasso_alpha;
        end
        params.result_codes = result_codes;
        params.kin_signals = kin_signals;
        params.train_start_idx = train_start_idx;
        params.train_end_idx = train_end_idx;
        cv_params = params; % these will be saved separately
        params.arrays = {'M1','PMd'};
        params.pert = pert;
        params.filedb = filedb;
        params.epochs = epochs;
        params.test_start_idx = test_start_idx;
        params.test_end_idx = test_end_idx;
        params.block_size_testing = block_size_testing;
        params.comment = params_comment;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load the trial data file
        filename = [filedb.Monkey{use_files(idx_file)} '_' filedb.Task{use_files(idx_file)} '_' filedb.Perturbation{use_files(idx_file)} '_' filedb.Date{use_files(idx_file)}];
        load(fullfile(rootDir,TDDir,[filename '.mat']));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
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
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % bin and prepare signals
        [trial_data, params] = glm_process_trial_data(trial_data,params);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Duplicating trial data for test set')
        trial_data_test = trimTD(trial_data,test_start_idx,test_end_idx);
        trial_data = trimTD(trial_data,test_start_idx,test_end_idx);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % filter trials and partition testing/training sets
        [~,trial_data] = getTDidx(trial_data,'result',result_codes);
        
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
        % Get potent/null space
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % project data into potent/null spaces
        for trial = 1:length(trial_data)
            temp = sqrt(trial_data(trial).PMd_spikes);
            if smooth_pca_proj_spikes
                temp = smoothSpikesForPCA(temp,params.bin_size,2*params.bin_size);
            end
            % de-mean data
            if center_pca
                temp = temp - repmat(params.pca_mu,size(temp,1),1);
            end
            temp = temp*params.pca_w;
            if ~ischar(params.pca_dims.PMd)
                temp = temp(:,params.pca_dims.PMd);
            end
            trial_data(trial).([cov_array '_potent']) = temp * params.V_potent;
            trial_data(trial).([cov_array '_null']) = temp * params.V_null;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % At last the data is processed!
        num_pred_vars = size(trial_data(1).PMd_potent,2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loop along neurons
        if do_cv
            disp('No CV file found or parameters did not match. Running CV...');
            [r2_full_cv,rr2_cv,r2_basic_cv] = deal(NaN(num_pred_vars,cv_folds,2));
        end
        bl_unit_models = repmat(struct('b_basic',[],'d_basic',[]),num_pred_vars,1);
        [r2_full,rr2,r2_basic] = deal(NaN(num_pred_vars,size(test_trials,1),2));
        
        x_basic = [ ...
            cat(1,trial_data(train_trials).pos), ...
            cat(1,trial_data(train_trials).vel), ...
            cat(1,trial_data(train_trials).speed), ...
            cat(1,trial_data(train_trials).force), ...
            cat(1,trial_data(train_trials).pos_rcb), ...
            cat(1,trial_data(train_trials).vel_rcb), ...
            cat(1,trial_data(train_trials).speed_rcb), ...
            cat(1,trial_data(train_trials).force_rcb)];
        x_full = [ ...
            x_basic, ...
            cat(1,trial_data(train_trials).PMd_null)];
        y_all = cat(1,trial_data(train_trials).PMd_potent);
        
        for var = 1:num_pred_vars
            tic;
            disp(['Running var ' num2str(var) ' of ' num2str(num_pred_vars)]);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Build baseline model
            y = y_all(:,var);
            
            % train basic model
            if ~isempty(x_basic)
            b_basic = [ones(size(y)) x_basic]\y;
            bl_unit_models(var).b_basic = b_basic;
            end
            b_full = [ones(size(y)) x_full]\y;
            bl_unit_models(var).b_full = b_full;
            
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
                    if ~isempty(x_basic)
                        b_basic = [ones(size(y_train)) x_basic(train_idx,:)]\y_train;
                    end
                    b_full = [ones(size(y_train)) x_full(train_idx,:)]\y_train;
                    
                    y_test = y(test_idx);
                    yfit_full = ([ones(length(test_idx),1), x_full(test_idx,:)]*b_full)';
                    if ~isempty(x_basic)
                        yfit_basic = ([ones(length(test_idx),1), x_basic(test_idx,:)]*b_basic)';
                    end
                    
                    % NEED TO UPDATE THIS
                    % compute R2
                    % get vector of bootstrap indices
                    bs = randi(length(y_test),length(y_test),num_bootstraps);
                    r2_full_cv(var,fold,:) = prctile(compute_vaf(y_test(bs),yfit_full(bs)),[2.5 97.5]);
                    if ~isempty(x_basic)
                        rr2_cv(var,fold,:) = prctile(compute_vaf(y_test(bs),yfit_basic(bs),yfit_full(bs)),[2.5 97.5]);
                        r2_basic_cv(var,fold,:) = prctile(compute_vaf(y_test(bs),yfit_basic(bs)),[2.5 97.5]);
                    end
                end
                disp(['R2 = ' num2str(mean(mean(r2_full_cv(var,:,:),3)))]);
                disp(['rR2 = ' num2str(mean(mean(rr2_cv(var,:,:),3)))]);
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Test throughout session on new data
            for e = 1:size(test_trials,1)
                
                x_basic = [ ...
                    cat(1,trial_data(test_trials(e,:)).pos), ...
                    cat(1,trial_data(test_trials(e,:)).vel), ...
                    cat(1,trial_data(test_trials(e,:)).speed), ...
                    cat(1,trial_data(test_trials(e,:)).force), ...
                    cat(1,trial_data(test_trials(e,:)).pos_rcb), ...
                    cat(1,trial_data(test_trials(e,:)).vel_rcb), ...
                    cat(1,trial_data(test_trials(e,:)).speed_rcb), ...
                    cat(1,trial_data(test_trials(e,:)).force_rcb)];
                x_full = [ ...
                    x_basic, ...
                    cat(1,trial_data(test_trials(e,:)).PMd_null)];
                y_all = cat(1,trial_data(test_trials(e,:)).PMd_potent);
                y = y_all(:,var);
                
                % fit models
                if ~isempty(x_basic)
                    yfit_basic = ([ones(size(x_basic,1),1), x_basic]*bl_unit_models(var).b_basic)';
                    trial_data_test(test_trials(e,:)).yfit_basic(:,var) = yfit_basic;
                end
                yfit_full = ([ones(size(x_full,1),1), x_full]*bl_unit_models(var).b_full)';
                trial_data_test(test_trials(e,:)).yfit_full(:,var) = yfit_full;
                
                
                % compute R2
                % get vector of bootstrap indices
                bs = randi(length(y),length(y),num_bootstraps);
                r2_full(var,e,:) = prctile(compute_vaf(y(bs),yfit_full(bs)),[2.5 97.5]);
                if ~isempty(x_basic)
                    rr2(var,e,:) = prctile(compute_vaf(y(bs),yfit_basic(bs),yfit_full(bs)),[2.5 97.5]);
                    r2_basic(var,e,:) = prctile(compute_vaf(y(bs),yfit_basic(bs)),[2.5 97.5]);
                end
                
            end
            toc;
        end % end neuron loop
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save results
        if do_cv
            save(fullfile(rootDir,resultsDir,outputSubdir,[pert '-null-potent_' filename '_cv.mat']),'r2_full_cv','r2_basic_cv','rr2_cv','cv_params');
        end
        
        results = struct( ...
            'r2_full', r2_full, ...
            'r2_basic', r2_basic, ...
            'rr2', rr2, ...
            'r2_full_cv', r2_full_cv, ...
            'r2_basic_cv', r2_basic_cv, ...
            'rr2_cv', rr2_cv, ...
            'bl_model', bl_unit_models);
        
        save(fullfile(rootDir,resultsDir,outputSubdir,[pert '-null-potent_' filename '_td.mat']),'trial_data_test','params');
        save(fullfile(rootDir,resultsDir,outputSubdir,[pert '-null-potent_' filename '.mat']),'results','params');
        clear file_results;
    end % end file loop
end % end perturbation loop
