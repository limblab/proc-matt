%%
% STUFF TO DO
%   - PCA as input
%   - null and potent space as input
%   - mahalonobis distance at time (e.g. target presentation, go cue, peak
%   speed, etc) between baseline and force as function of adaptation (for
%   PCA of M1, PMd, etc)
%
%
clear;
clc;
close all;

dataSummary;

outputSubdir = 'all';
params_file = ''; % optionally give .mat file containing params to use
% Must be in the outputSubdir

sessions = { ...
%     'Chewie','2016-09-15'; ...
    'Chewie','2016-09-19'; ...
    'Chewie','2016-10-05'; ...
    'Chewie','2016-10-07'; ...
%     'Mihili','2014-02-17'; ...
%     'Mihili','2014-02-18'; ...
%     'Mihili','2014-03-07'; ...
%     'Mihili','2015-06-15'; ...
%     'Mihili','2015-06-16'; ...
%     'Mihili','2015-06-17'; ...
    };

monkeys = unique(sessions(:,1));
tasks = {'CO'};
perts = {'FF'};
dates = sessions(:,2);
% {COVARIATE, PREDICTED}, will loop along rows
% array_pairs = {'PMd','M1';'PMd','PMd';'M1','M1'};
array_pairs = {'PMd','PMd'};

%%
if isempty(params_file)
    dt                  = 0.01; % time step size for data
    bin_size            = 5;    % how many samples to group together when rebinning (bin width is num_samples*dt)
    
    result_codes        = {'R'}; % which trials to include
    
    training_data       = {'BL','all'};
    testing_data        = {'AD','all'; ...
                           'WO','all'};
    block_size_testing  = 1; % size of blocks in # trials for AD/WO
    % note: can change this and still reload CV
    
    % how to truncate trials {idx name, number of bins after}
    train_start_idx     = {'idx_target_on',0}; %{'idx_target_on',-4}
    train_end_idx       = {'idx_trial_end',-2}; %{'idx_go_cue',4}
    test_start_idx      = {'idx_target_on',0}; % NOTE: CAN CHANGE THESE
    test_end_idx        = {'idx_trial_end',-2}; %   AND STILL RELOAD CV
    %   NOTE: this is after rebinning at the moment
    
    cv_folds            = 10; % how many folds for cross-validation
    cv_rand             = true; % whether to randomly select CV datapoints
    
    block_size_fr_test  = 100;   % how many trials to group for FR test
    fr_test_alpha       = 1e-3; % p value cut off for t-test
    fr_min              = 0.3; % minimum session-wide spiking for inclusion
    
    do_lasso            = false;  % if you want to regularize
    lasso_lambda        = 0.0083; % lambda parameter for lasso
    lasso_alpha         = 0.01;   % alpha parameter for lasso
    
    do_pca              = false; % include PCA as covariates instead of spikes
    pca_dims            = 'all'; % 'all' or specify which dimensions, e.g. 1:30
    
    do_rcb              = true;  % use raised cosine basis functions for history
    do_all_history      = false; % include history for all units
    do_self_history     = false;
    unit_lags           = 2;     % how many bins in the past for each neuron to include as covariates (if ~do_rcb)
    rcb_hpeaks          = [dt*bin_size,0.1];
    rcb_b               = 0.1;
    %   NOTE: this is after rebinning at the moment
    
    kin_signals         = {'pos','vel','speed'}; % names of kinematic covariates
    kin_lags            = 2; % how many bins to lag (if ~do_rcb) or how many raised cosines (if do_rcb)
    %   NOTE: this is after rebinning at the moment
    
    num_bootstraps      = 1000; % how many resamples for R2
    
else % load up params from whatever file
    load(params_file,'params');
    fn = fieldnames(params);
    for i = 1:length(fn)
        eval([fn{i} ' = params.' fn{i} ';']);
    end, clear i good_cells cov_array pred_array pert epochs pca_w;
end

%% Here is where I check all of the inputs
if bin_size < 1 || block_size_testing < 1
    error('Invalid bin size or testing block size');
elseif bin_size == 1
    disp('---------------------------------------------------------')
    disp('BIN SIZE IS 1... SOMETIMES THIS CRASHES THINGS. CONTINUE?');
    disp('---------------------------------------------------------')
    pause;
end

%%
if ~exist(fullfile(rootDir,TDDir,outputSubdir),'dir')
    mkdir(fullfile(rootDir,TDDir,outputSubdir));
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
            if do_pca
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
            params.train_start_idx = train_start_idx;
            params.train_end_idx = train_end_idx;
            cv_params = params; % these will be saved separately
            params.pert = pert;
            params.filedb = filedb;
            params.epochs = epochs;
            params.test_start_idx = test_start_idx;
            params.test_end_idx = test_end_idx;
            params.block_size_testing = block_size_testing;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load the trial data file
            filename = [filedb.Monkey{use_files(idx_file)} '_' filedb.Task{use_files(idx_file)} '_' filedb.Perturbation{use_files(idx_file)} '_' filedb.Date{use_files(idx_file)}];
            load(fullfile(rootDir,TDDir,[filename '.mat']));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % will save cross validated BL predictions so I don't have to
            % redo them to to apply the model to other time points
            cv_fn = fullfile(rootDir,TDDir,outputSubdir,[pert '-' cov_array '-' pred_array '_' filename '_cv.mat']);
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
            % Calculate speed if desired
            if any(ismember(kin_signals,'speed'))
                for trial = 1:length(trial_data)
                    trial_data(trial).speed = hypot(trial_data(trial).vel(:,1), ...
                        trial_data(trial).vel(:,2));
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % filter out neurons
            trial_data = getCommonUnits(trial_data);
            
            arrays = unique(array_pairs(idx_arrays,:));
            good_cells = cell(1,length(arrays));
            for array = 1:length(arrays);
                % make sure firing is significantly non-zero across the
                % whole session. Usually this means there is a sort problem
                % or there is noise and I lose a neuron, a real neuron
                % wouldn't completely shut off.
                all_fr=cell2mat(cellfun(@(x) sqrt(sum(x,1))', ...
                    cellfun(@(x) full(x),{trial_data.([arrays{array} '_spikes'])}, ...
                    'Uni',0),'Uni',0));
                
                % additionally, make sure the baseline firing (with the
                % model fit) is sufficiently high. I do this separately
                % from the last one because I want to allow neurons to
                % increase or decrease firing as they see fit during
                % learning
                bl_fr = [];
                trial_idx = find(strcmpi({trial_data.epoch},'BL'));
                for iTrial = trial_idx
                    idx = trial_data(iTrial).(train_start_idx{1})+train_start_idx{2}:trial_data(iTrial).(train_end_idx{1})+train_end_idx{2};
                    temp = trial_data(iTrial).([arrays{array} '_spikes']);
                    bl_fr = [bl_fr, (sqrt(sum(temp(idx,:),1))./(size(temp,1)*dt))'];
                end
                
                % ensure all blocks are significantly non-zero
                blocks = 1:block_size_fr_test:size(all_fr,2);
                p = zeros(size(all_fr,1),length(blocks));
                for block = 2:length(blocks)
                    for unit = 1:size(all_fr,1)
                        [~,p(unit,block)] = ttest(all_fr(unit,blocks(block-1):blocks(block)),0,'tail','right');
                    end
                end
                good_cells{array} = find(all(p <= fr_test_alpha,2) & mean(bl_fr,2) >= fr_min)';
                
                disp(['Removing ' num2str(size(all_fr,1) - length(good_cells{array})) ' low-firing cells...']);
                
                % take bad cells out of trial_data
                for trial = 1:length(trial_data)
                    temp = trial_data(trial).([arrays{array} '_spikes']);
                    trial_data(trial).([arrays{array} '_spikes']) = temp(:,good_cells{array});
                    temp = trial_data(trial).([arrays{array} '_unit_guide']);
                    trial_data(trial).([arrays{array} '_unit_guide']) = temp(good_cells{array},:);
                end
            end, clear all_fr bl_fr num_blocks i array r s temp blocks;
            params.good_cells = good_cells;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Re-bin data
            if bin_size > 1
                disp('Binning...');
                trial_data = truncateAndBin(trial_data,bin_size);
                params.dt = dt*bin_size;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if do_rcb
                % Convolve with basis functions before re-binning
                disp('Convolving with raised cosines...');
                rcb_which_vars = cell(1,length(arrays));
                for i = 1:length(arrays)
                    rcb_which_vars{i} = [arrays{i} '_spikes'];
                end
                rcb_which_vars = [rcb_which_vars, kin_signals];
                params.rcb_which_vars = rcb_which_vars;
                tic;
                trial_data = convBasisFunc(trial_data,params);
                toc;
                
                for i = 1:length(kin_signals)
                    % shift kinematics by kin_lags backwards
                    for j = 1:length(trial_data)
                        temp = trial_data(j).(kin_signals{i});
                        trial_data(j).(kin_signals{i}) = [temp(kin_lags+1:end,:); NaN(kin_lags,size(temp,2))];
                        temp = trial_data(j).([kin_signals{i} '_shift']);
                        trial_data(j).([kin_signals{i} '_shift']) = [temp(kin_lags+1:end,:); NaN(kin_lags,size(temp,2))];
                    end
                end                
            else
                error('FIX KINEMATIC SHIFT');
                if unit_lags > 0
                    % Duplicate and shift
                    build_inputs = cell(1,2*(length(kin_signals) + length(arrays)));
                    for i = 1:length(arrays)
                        build_inputs{(i-1)*2+1} = [arrays{i} '_spikes'];
                        build_inputs{(i-1)*2+2} = unit_lags;
                    end
                    if ~isempty(kin_signals)
                        for i = 1:length(kin_signals)
                            build_inputs{2*length(arrays) + (i-1)*2+1} = kin_signals{i};
                            build_inputs{2*length(arrays) + (i-1)*2+2} = max(kin_lags);
                        end
                    end
                    trial_data = dupeAndShift(trial_data,build_inputs);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % filter trials and partition testing/training sets
            if isfield(trial_data,'result') % old format didn't include result and only had 'R'
                disp('Result not found. All trials are reward trials.');
                trial_data = trial_data(ismember({trial_data.result},result_codes));
            end
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
            if ~strcmpi(cov_array,pred_array) && do_pca % we don't have to leave one out so just fit once
                w = getPCA(trial_data,struct('array',cov_array,'bin_size',dt*bin_size,'trial_idx',1:length(trial_data),'neurons',1:num_cov_neurons));
                params.pca_w = w;
                % do SVD to get null/potent spaces
                
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
                
                if strcmpi(cov_array,pred_array) && do_pca  % recompute PCA each time
                    % leave predicted neuron out of PCA
                    neuron_idx = 1:num_pred_neurons; neuron_idx(unit) = [];
                    w = getPCA(trial_data,struct('array',cov_array,'bin_size',dt*bin_size,'trial_idx',1:length(trial_data),'neurons',neuron_idx));
                    % used in glm_prep_inputs
                    params.pca_w = w;
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
                    
                    % compute R2
                    % get vector of bootstrap indices
                    bs = randi(length(y),length(y),num_bootstraps);
                    pr2_full(unit,e,:) = prctile(compute_pseudo_R2(y(bs),yfit_full(bs),mean(y)),[2.5 97.5]);
                    if ~isempty(b_basic)
                        rpr2(unit,e,:) = prctile(compute_rel_pseudo_R2(y(bs),yfit_basic(bs),yfit_full(bs)),[2.5 97.5]);
                        pr2_basic(unit,e,:) = prctile(compute_pseudo_R2(y(bs),yfit_basic(bs),mean(y)),[2.5 97.5]);
                    end
                    
                    %                     rpr2(unit,e,:) = bootci(num_bootstraps,{@compute_rel_pseudo_R2,y, yfit_basic', yfit_full'});
                    %                     pr2_full(unit,e,:) = bootci(num_bootstraps,{@compute_pseudo_R2,y, yfit_full', mean(pred_fr(:,i))});
                    %                     pr2_basic(unit,e,:) = bootci(num_bootstraps,{@compute_pseudo_R2,y, yfit_basic', mean(pred_fr(:,i))});
                    
                    if 0
                        figure;
                        hold all;
                        plot(y);
                        plot(yfit_full);
                        plot(yfit_basic);
                        legend({['rpr2=' num2str(mean(rpr2(unit,e,:)))],['pr2_f=' num2str(mean(pr2_full(unit,e,:)))],['pr2_b=' num2str(mean(pr2_basic(unit,e,:)))]});
                        pause;
                        close all;
                    end
                end
                toc;
%                 subplot(1,2,1); hold all; plot(mean(pr2_full(unit,:,:),3));
%                 subplot(1,2,2); hold all; plot(mean(rpr2(unit,:,:),3));
%                 drawnow;
            end % end neuron loop
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Save results
            if do_cv
                save(fullfile(rootDir,TDDir,outputSubdir,[pert '-' cov_array '-' pred_array '_' filename '_cv.mat']),'pr2_full_cv','pr2_basic_cv','rpr2_cv','cv_params');
            end
            
            results = struct( ...
                'pr2_full', pr2_full, ...
                'pr2_basic', pr2_basic, ...
                'rpr2', rpr2, ...
                'pr2_full_cv', pr2_full_cv, ...
                'pr2_basic_cv', pr2_basic_cv, ...
                'rpr2_cv', rpr2_cv, ...
                'bl_model', bl_unit_models);
            
            save(fullfile(rootDir,TDDir,outputSubdir,[pert '-' cov_array '-' pred_array '_' filename '.mat']),'results','params');
            clear file_results;
        end % end array pair loop
    end % end file loop
end % end perturbation loop
