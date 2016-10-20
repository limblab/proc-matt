%% THIS ONE WILL BUILD ON ADAPTATION WITH FORCE AND TEST ON BASELINE
clear;
clc;
close all;

dataSummary;

outputSubdir = '50msec_bins_test';

monkeys = {'Chewie'};
tasks = {'CO'};
perts = {'FF'};
dates = {'2016-10-07'};

result_codes = {'R'};

% {COVARIATE, PREDICTED}, will loop along rows
% array_pairs = {'M1','M1'; ...
% 'PMd','PMd'; ...
array_pairs = {'PMd','M1'; 'PMd','PMd';'M1','M1'};

% how to break up the epochs into blocks
num_blocks_bl = 0; % how many BL blocks to exclude for testing
block_size_testing = 20; % size of blocks in # trials

dt = 0.01; % time step size for data
num_samples = 5; % how many samples to group together when rebinning (bin width is num_samples*dt)

block_size_fr_test = 20;
fr_test_alpha = 1e-4; % p value cut off for t-test
fr_min = 2; % minimum session-wide spiking for inclusion

% {idx name, number of bins after}
start_idx = {'idx_target_on',-1};
end_idx = {'idx_trial_end',-2};
%   NOTE: this is after rebinning at the moment

do_lasso = false;
lasso_lambda = 0.0083;
lasso_alpha = 0.01;

do_pca = false;

do_all_history = false;
unit_lags = 3;% how many bins in the past for each neuron to include as covariates (duplicate and shift)
%   NOTE: this is after rebinning at the moment

do_cv = true;
num_folds = 10; % how many folds for cross-validation

do_kin_model = false;
kin_signals = {'pos','vel'};
kin_lags = 2; % how many bins to lag
%   NOTE: this is after rebinning at the moment

bootstrap_r2 = true; % if model fit is cross validated this isn't as necessary because I can use that for "significance"
num_bootstraps = 1000;

% for bootstrapping
opts.UseParallel = true;


%%
if ~exist(fullfile(rootDir,TDDir,outputSubdir),'dir')
    mkdir(fullfile(rootDir,TDDir,outputSubdir));
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% disp(['Training on ' training_epoch{1} '...']);

for idx_pert = 1:length(perts)
    pert = perts{idx_pert};
    disp(['Starting perturbation ' pert '...']);
    
    for idx_arrays = 1:size(array_pairs,1)
        
        disp(['STARTING ARRAY PAIR ' array_pairs{idx_arrays,1} ' -> ' array_pairs{idx_arrays,2} '...']);
        
        cov_array = array_pairs{idx_arrays,1};
        pred_array = array_pairs{idx_arrays,2};
        
        use_date_idx = ismember(filedb.Monkey,monkeys) & ismember(filedb.Perturbation,perts(idx_pert)) & ismember(filedb.Task,tasks);
        if ~isempty(dates)
            use_date_idx = use_date_idx & ismember(filedb.Date,dates);
        end
        use_files = find(use_date_idx);
        
        clear results;
        for idx_file = 1:length(use_files)
            disp(['File ' num2str(idx_file) ' of ' num2str(length(use_files)) '...'])
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Package up parameters for saving
            clear params;
            params.pert = pert;
            params.result_codes = result_codes;
            params.filedb = filedb;
            params.use_files = use_files;
            params.cov_array = cov_array;
            params.pred_array = pred_array;
            params.num_blocks_bl = num_blocks_bl;
            params.block_size_testing = block_size_testing;
            params.block_size_fr_test = block_size_fr_test;
            params.fr_test_alpha = fr_test_alpha;
            params.fr_min = fr_min;
            params.dt = dt;
            params.num_samples = num_samples;
            params.num_folds = num_folds;
            params.unit_lags = unit_lags;
            params.do_lasso = do_lasso;
            params.do_all_history = do_all_history;
            params.do_kin_model = do_kin_model;
            params.bootstrap_r2 = bootstrap_r2;
            params.num_bootstraps = num_bootstraps;
            params.lasso_lambda = lasso_lambda;
            params.lasso_alpha = lasso_alpha;
            params.start_idx = start_idx;
            params.end_idx = end_idx;
            params.kin_signals = kin_signals;
            params.kin_lags = kin_lags;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Load the trial data file
            filename = [filedb.Monkey{use_files(idx_file)} '_' filedb.Task{use_files(idx_file)} '_' filedb.Perturbation{use_files(idx_file)} '_' filedb.Date{use_files(idx_file)}];
            load(fullfile(rootDir,TDDir,[filename '.mat']));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Re-bin data
            if num_samples > 1
                trial_data = truncateAndBin(trial_data,num_samples);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % filter out neurons
            arrays = unique(array_pairs(idx_arrays,:));
            good_cells = cell(1,length(arrays));
            for j = 1:length(arrays);
                % make sure firing is significantly non-zero across the
                % whole session. Usually this means there is a sort problem
                % or there is noise and I lose a neuron, a real neuron
                % wouldn't completely shut off.
                all_fr=cell2mat(cellfun(@(x) sqrt(sum(x,1))', ...
                    cellfun(@(x) full(x),{trial_data.([arrays{j} '_spikes'])}, ...
                    'UniformOutput',false),'UniformOutput',false));
                
                % additionally, make sure the baseline firing (with the
                % model fit) is sufficiently high. I do this separately
                % from the last one because I want to allow neurons to
                % increase or decrease firing as they see fit during
                % learning
                bl_fr = cell2mat(cellfun(@(x) sqrt(sum(x,1))', ...
                    cellfun(@(x) full(x),{trial_data(strcmpi({trial_data.epoch},'bl')).([arrays{j} '_spikes'])}, ...
                    'UniformOutput',false),'UniformOutput',false));
                
                % ensure all blocks are significantly non-zero
                blocks = 1:block_size_fr_test:size(all_fr,2);
                p = zeros(size(all_fr,1),length(blocks));
                for k = 2:length(blocks)
                    for i = 1:size(all_fr,1)
                        [~,p(i,k)] = ttest(all_fr(i,blocks(k-1):blocks(k)),0,'tail','right');
                    end
                end
                good_cells{j} = find(all(p <= fr_test_alpha,2) & mean(bl_fr,2) >= fr_min)';
                
                disp(['Removing ' num2str(size(all_fr,1) - length(good_cells{j})) ' low-firing cells...']);
                
                % take bad cells out of trial_data
                for i = 1:length(trial_data)
                    temp = trial_data(i).([arrays{j} '_spikes']);
                    trial_data(i).([arrays{j} '_spikes']) = temp(:,good_cells{j});
                    temp = trial_data(i).([arrays{j} '_unit_guide']);
                    trial_data(i).([arrays{j} '_unit_guide']) = temp(good_cells{j},:);
                end
            end, clear all_fr bl_fr num_blocks i j r s temp arrays blocks;
            params.good_cells = good_cells;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Duplicate and shift
            trial_data = getCommonUnits(trial_data);
            trial_data = dupeAndShift(trial_data,'pos',max(kin_lags),'vel',max(kin_lags),'M1_spikes',unit_lags,'PMd_spikes',unit_lags);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Truncate to trial markers
            trial_data = truncateAndBin(trial_data,start_idx,end_idx);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % filter out trials
            trial_data = trial_data(ismember({trial_data.result},result_codes));
            
            training_idx = find(strcmpi({trial_data.epoch},'BL'));
            testing_idx_bl = training_idx(end-num_blocks_bl*block_size_testing+1:end);
            training_idx = training_idx(1:end-num_blocks_bl*block_size_testing);
            
            % get testing indices throughout AD and WO
            testing_idx_ad = (1:sum(strcmpi({trial_data.epoch},'AD')))';
            testing_idx_ad = reshape(testing_idx_ad(1:end-rem(length(testing_idx_ad),block_size_testing)),block_size_testing,(length(testing_idx_ad) - rem(length(testing_idx_ad),block_size_testing))/block_size_testing)';
            testing_idx_wo = (1:sum(strcmpi({trial_data.epoch},'WO')))';
            testing_idx_wo = reshape(testing_idx_wo(1:end-rem(length(testing_idx_wo),block_size_testing)),block_size_testing,(length(testing_idx_wo) - rem(length(testing_idx_wo),block_size_testing))/block_size_testing)';
            
            bl_inds = find(strcmpi({trial_data.epoch},'bl'));
            ad_inds = find(strcmpi({trial_data.epoch},'ad'));
            wo_inds = find(strcmpi({trial_data.epoch},'wo'));
            
            testing_idx = [bl_inds(testing_idx_bl); ...
                           ad_inds(testing_idx_ad); ...
                           wo_inds(testing_idx_wo)];
            clear testing_idx_bl testing_idx_ad testing_idx_wo bl_inds ad_inds wo_inds;
            
            params.training_idx = training_idx;
            params.testing_idx = testing_idx;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get training covariates
            [cov_fr, cov_fr_shift, pred_fr, pred_fr_shift, cov_kin, cov_kin_shift] = get_covariates(trial_data,training_idx,params);
            %%% get full Baseline model
            disp('Fitting training model...');
            bl_unit_models = repmat(struct('b_basic',[],'b_full',[]),size(pred_fr,2),1);
            tic;
            for i = 1:size(pred_fr,2) % loop along neurons
                [x_self, x_kin] = deal([]);
                
                idx = true(1,size(cov_fr,2));
                if strcmpi(cov_array,pred_array),
                    % take out this neuron and compute the new PCA
                    idx(i) = false;
                    
                end
                
                y = pred_fr(:,i);
                
                if do_pca % get PCA space
                    % get for covariate array
                    %[trial_data,w,scores,eigen] = getPCA(trial_data, struct('array',cov_array,  'bin_size',dt*num_samples, 'neurons',find(idx)));
                else
                    x_fr = cov_fr(:,idx);
                end
                
                % add shifted data as covariates if desired
                if unit_lags > 0
                    % add time shifted other neurons
                    if do_all_history && ~do_pca
                        for j = 1:unit_lags, x_fr = [x_fr, cov_fr_shift(:,idx+(j-1)*size(cov_fr,2))]; end
                    end
                    x_self = zeros(size(pred_fr,1),unit_lags);
                    for j = 1:unit_lags, x_self(:,j) = pred_fr_shift(:,i+(j-1)*size(pred_fr,2)); end
                end
                
                if ~isempty(kin_signals)
                    for j = 1:length(kin_lags)
                        if kin_lags(j) == 0
                            x_kin = [x_kin, cov_kin];
                        else
                            x_kin = [x_kin, cov_kin_shift(:,1+(kin_lags(j)-1)*size(cov_kin,2):kin_lags(j)*size(cov_kin,2))];
                        end
                    end
                end
                
                x_full = [x_self, x_kin, x_fr];
                x_basic = [x_self, x_kin];
                
                if do_lasso
                    [b,s] = lassoglm(x_full,y,'poisson','lambda',lasso_lambda,'alpha',lasso_alpha);
                    b_full = [s.Intercept; b];
                    [b,s] = lassoglm(x_basic,y,'poisson','lambda',lasso_lambda,'alpha',lasso_alpha);
                    b_basic = [s.Intercept; b];
                else
                    b_full = glmfit(x_full,y,'poisson');
                    b_basic = glmfit(x_basic,y,'poisson');
                end
                
                bl_unit_models(i).b_basic = b_basic;
                bl_unit_models(i).b_full = b_full;
            end, clear x x_self x_kin x_fr b_basic b_full y j idx;
            
            if do_kin_model,
                error('Kin model currently broken.');
                x_fr = cov_fr;
                % add shifted data as covariates if desired
                if unit_lags > 0 && do_all_history
                    % add time shifted other neurons
                    for j = 1:unit_lags, x_fr = [x_fr, squeeze(cov_fr_shift(:,:,j))]; end
                end
                
                bl_b_kin = zeros(size(x_fr,2)+1,size(cov_kin,2));
                for i = 1:size(cov_kin,2)
                    if do_lasso
                        [b,s] = lassoglm(x_fr,cov_kin(:,i),'normal');
                        bl_b_kin(:,i) = [s.Intercept; b];
                    else
                        bl_b_kin(:,i) = glmfit(x_fr,cov_kin(:,i),'normal');
                    end
                end
            else
                bl_b_kin = [];
            end
            
            toc;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% run cross-validated model on baseline data to assess performance
            if do_cv
                disp('Cross-validating training model...');
                [pr2_full_cv,rpr2_cv,pr2_basic_cv] = deal(zeros(size(pred_fr,2),num_folds,2));
                
                tic;
                for i = 1:size(pred_fr,2)
                    disp(['Neuron ' num2str(i) ' of ' num2str(size(pred_fr,2))]);
                    
                    [x_self, x_kin] = deal([]);
                    
                    idx = true(1,size(cov_fr,2));
                    if strcmpi(cov_array,pred_array), idx(i) = false; end
                    
                    x_fr = cov_fr(:,idx);
                    % add shifted data as covariates if desired
                    if unit_lags > 0
                        if do_all_history
                            % add time shifted other neurons
                            for j = 1:unit_lags, x_fr = [x_fr, cov_fr_shift(:,idx+(j-1)*size(cov_fr,2))]; end
                        end
                        % now add self-history
                        x_self = zeros(size(pred_fr,1),unit_lags);
                        for j = 1:unit_lags, x_self(:,j) = pred_fr_shift(:,i+(j-1)*size(pred_fr,2)); end
                    end
                    
                    if ~isempty(kin_signals)
                        for j = 1:length(kin_lags)
                            if kin_lags(j) == 0
                                x_kin = [x_kin, cov_kin];
                            else
                                x_kin = [x_kin, cov_kin_shift(:,1+(kin_lags(j)-1)*size(cov_kin,2):kin_lags(j)*size(cov_kin,2))];
                            end
                        end
                    end
                    
                    x_full = [x_self, x_kin, x_fr];
                    x_basic = [x_self, x_kin];
                    
                    folds = linspace(1,size(cov_fr,1),num_folds+1);
                    for j = 1:num_folds
                        test_idx = ceil(folds(j)):floor(folds(j+1));
                        train_idx = true(1,size(cov_fr,1)); train_idx(test_idx) = false;
                        
                        y = pred_fr(train_idx,i);
                        if do_lasso
                            [b,s] = lassoglm(x_full(train_idx,:),y,'poisson','lambda',lasso_lambda,'alpha',lasso_alpha);
                            b_full = [s.Intercept; b];
                            [b,s] = lassoglm(x_basic(train_idx,:),y,'poisson','lambda',lasso_lambda,'alpha',lasso_alpha);
                            b_basic = [s.Intercept; b];
                        else
                            b_full = glmfit(x_full(train_idx,:),y,'poisson');
                            b_basic = glmfit(x_basic(train_idx,:),y,'poisson');
                        end
                        
                        y = pred_fr(test_idx,i);
                        yfit_full = exp([ones(length(test_idx),1), x_full(test_idx,:)]*b_full)';
                        yfit_basic = exp([ones(length(test_idx),1), x_basic(test_idx,:)]*b_basic)';
                        
                        % compute R2
                        try
                            rpr2_cv(i,j,:) = bootci(num_bootstraps,{@compute_rel_pseudo_R2,y, yfit_basic', yfit_full'},'Options',opts);
                            pr2_full_cv(i,j,:) = bootci(num_bootstraps,{@compute_pseudo_R2,y, yfit_full', mean(pred_fr(train_idx,i))},'Options',opts);
                            pr2_basic_cv(i,j,:) = bootci(num_bootstraps,{@compute_pseudo_R2,y, yfit_basic', mean(pred_fr(train_idx,i))},'Options',opts);
                        catch
                            rpr2_cv(i,j,:) = nan(1,2);
                            pr2_full_cv(i,j,:) = nan(1,2);
                            pr2_basic_cv(i,j,:) = nan(1,2);
                        end
                    end
                end, clear i j y x_full x_basic b_full b_basic yfit_full yfit_basic y train_idx test_idx;
                
                %%%%%%%%%%%%%%%% DO kinematic predictions
                if do_kin_model
                    error('Kin model is currently broken.');
                    pr2_kin_cv = zeros(size(cov_kin,2),num_folds,2);
                    x_fr = cov_fr;
                    % add shifted data as covariates if desired
                    if unit_lags > 0 && do_all_history
                        % add time shifted other neurons
                        for j = 1:unit_lags, x_fr = [x_fr, squeeze(cov_fr_shift(:,:,j))]; end
                    end
                    
                    for i = 1:size(cov_kin,2)
                        folds = linspace(1,size(cov_fr,1),num_folds+1);
                        for j = 1:num_folds
                            test_idx = ceil(folds(j)):floor(folds(j+1));
                            train_idx = true(1,size(cov_fr,1)); train_idx(test_idx) = false;
                            
                            if do_lasso
                                [b,s] = lassoglm(x_fr(train_idx,:),cov_kin(train_idx,i),'normal');
                                b_kin(:,k) = [s.Intercept; b];
                            else
                                b_kin = glmfit(x_fr(train_idx,:),cov_kin(train_idx,i),'normal');
                            end
                            
                            y = cov_kin(test_idx,i);
                            yfit_kin = ([ones(length(test_idx),1), x_fr(test_idx,:)]*b_kin)';
                            
                            try
                                pr2_kin_cv(i,j,:) = bootci(num_bootstraps,{@CalculateR2,y, yfit_kin'},'Options',opts);
                            catch
                                pr2_kin_cv(i,j,:) = nan(1,2);
                            end
                            
                        end
                    end
                else
                    pr2_kin_cv = [];
                end % end do_kin_model
                toc;
            else
                rpr2_cv = [];
                pr2_full_cv = [];
                pr2_basic_cv = [];
                pr2_kin_cv = [];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% now use training model to predict spiking in each epoch
            disp('Calculating metrics...');
            
            [pr2_full,rpr2,pr2_basic] = deal(zeros(size(pred_fr,2),size(testing_idx,1),2));
            if do_kin_model, pr2_kin = zeros(size(bl_b_kin,2),size(testing_idx,1),2); end
            
            % loop along the blocks of test trials
            for e = 1:size(testing_idx,1)
                disp(['Predicting ' num2str(e) ' of ' num2str(size(testing_idx,1)) '...']);
                
                [cov_fr, cov_fr_shift, pred_fr, pred_fr_shift, cov_kin, cov_kin_shift] = get_covariates(trial_data, testing_idx(e,:), params);
                
                tic;
                for i = 1:size(pred_fr,2)
                    [x_self, x_kin] = deal([]);
                    
                    b_basic = bl_unit_models(i).b_basic;
                    b_full = bl_unit_models(i).b_full;
                    
                    idx = true(1,size(cov_fr,2));
                    if strcmpi(cov_array,pred_array), idx(i) = false; end
                    
                    y = pred_fr(:,i);
                    
                    x_fr = cov_fr(:,idx);
                    % add shifted data as covariates if desired
                    if unit_lags > 0
                        if do_all_history % add time shifted other neurons
                            for j = 1:unit_lags, x_fr = [x_fr, cov_fr_shift(:,idx+(j-1)*size(cov_fr,2))]; end
                        end
                        
                        % now add self-history
                        x_self = zeros(size(pred_fr,1),unit_lags);
                        for j = 1:unit_lags, x_self(:,j) = pred_fr_shift(:,i+(j-1)*size(pred_fr,2)); end
                    end
                    
                    if ~isempty(kin_signals)
                        for j = 1:length(kin_lags)
                            if kin_lags(j) == 0
                                x_kin = [x_kin, cov_kin];
                            else
                                x_kin = [x_kin, cov_kin_shift(:,1+(kin_lags(j)-1)*size(cov_kin,2):kin_lags(j)*size(cov_kin,2))];
                            end
                        end
                    end
                    
                    x_full = [x_self, x_kin, x_fr];
                    x_basic = [x_self, x_kin];
                    
                    yfit_full = exp([ones(size(x_full,1),1), x_full]*b_full)';
                    yfit_basic = exp([ones(size(x_basic,1),1), x_basic]*b_basic)';                  
                    
                    % compute R2
                    try
                        rpr2(i,e,:) = bootci(num_bootstraps,{@compute_rel_pseudo_R2,y, yfit_basic', yfit_full'},'Options',opts);
                        pr2_full(i,e,:) = bootci(num_bootstraps,{@compute_pseudo_R2,y, yfit_full', mean(pred_fr(:,i))},'Options',opts);
                        pr2_basic(i,e,:) = bootci(num_bootstraps,{@compute_pseudo_R2,y, yfit_basic', mean(pred_fr(:,i))},'Options',opts);
                    catch
                        rpr2(i,e,:) = nan(1,2);
                        pr2_full(i,e,:) = nan(1,2);
                        pr2_basic(i,e,:) = nan(1,2);
                    end
                end
                
                %%%%%%%% DO KINEMATIC PREDICTIONS
                if do_kin_model
                    error('Kin model is currently broken.');
                    x_fr = cov_fr;
                    % add shifted data as covariates if desired
                    if unit_lags > 0 && do_all_history % add time shifted other neurons
                        for j = 1:unit_lags, x_fr = [x_fr, squeeze(cov_fr_shift(:,idx,j))]; end
                    end
                    
                    for i = 1:size(bl_b_kin,2)
                        b_kin = bl_b_kin(:,i);
                        y = cov_kin(:,i);
                        yfit_kin = ([ones(size(x_fr,1),1), x_fr]*b_kin)';
                        
                        try
                            pr2_kin(i,e,:) = bootci(num_bootstraps,{@CalculateR2,y, yfit_kin'},'Options',opts);
                        catch
                            pr2_kin(i,e,:) = nan(1,2);
                        end
                    end
                else
                    pr2_kin = [];
                end
                
                toc;
            end

            results(idx_file).pr2_full = pr2_full;
            results(idx_file).pr2_basic = pr2_basic;
            results(idx_file).rpr2 = rpr2;
            results(idx_file).pr2_full_cv = pr2_full_cv;
            results(idx_file).pr2_basic_cv = pr2_basic_cv;
            results(idx_file).rpr2_cv = rpr2_cv;
            results(idx_file).pr2_kin = pr2_kin;
            results(idx_file).pr2_kin_cv = pr2_kin_cv;
            results(idx_file).bl_model = bl_unit_models;
            results(idx_file).bl_b_kin = bl_b_kin;
            
            
            file_results.pr2_full = pr2_full;
            file_results.pr2_basic = pr2_basic;
            file_results.rpr2 = rpr2;
            file_results.pr2_full_cv = pr2_full_cv;
            file_results.pr2_basic_cv = pr2_basic_cv;
            file_results.rpr2_cv = rpr2_cv;
            file_results.pr2_kin = pr2_kin;
            file_results.pr2_kin_cv = pr2_kin_cv;
            file_results.bl_model = bl_unit_models;
            file_results.bl_b_kin = bl_b_kin;
            save(fullfile(rootDir,TDDir,outputSubdir,[pert '-' cov_array '-' pred_array '_' filename '.mat']),'file_results','params');
            clear file_results;
        end
        
        save(fullfile(rootDir,TDDir,outputSubdir,[pert '-' cov_array '-' pred_array '.mat']),'results','params');
    end
end, clear use_date_idx sessionList idx_arrays idx_pert;
