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
    
    load(fullfile(rootDir,TDDir,'glm_encoding_pca_without_smoothing','trainad_null_nokin',['FF-PMd-M1_' filename '_cv.mat']),'cv_params');
    params = cv_params;
    params.arrays = {'M1','PMd'};
    params.do_rcb = false;
    [trial_data, params] = glm_process_trial_data(trial_data,params);
    load(fullfile(rootDir,TDDir,'glm_encoding_pca_without_smoothing','trainad_null_nokin',['FF-PMd-M1_' filename '.mat']),'params');
    
    if isfield(trial_data,'result') % old format didn't include result and only had 'R'
        trial_data = trial_data(ismember({trial_data.result},result_codes));
    else
        disp('Result not found. All trials are reward trials.');
    end
    
    pred_array = 'M1';
    cov_array = 'PMd';
    params.pred_array = pred_array;
    params.cov_array = cov_array;
    
    if ischar(neuron_id)
        neuron_id = trial_data(1).([pred_array '_unit_guide']);
    end
    
    % prep trial_data for new fields
    for trial = 1:length(trial_data)
        trial_data(trial).Pred_spikes = zeros(size(trial_data(trial).([pred_array '_spikes'])));
        %         trial_data(trial).pred_spikes_pr2 = zeros(size(trial_data(trial).([pred_array '_spikes']),2),2);
        %     trial_data(trial).pred_spikes_rpr2 = zeros(size(trial_data(trial).([pred_array '_spikes']),2),2);
    end
    
    % loop along trials and neurons and build predicted spikes
    if 0
        pred_models = {'potent','null'};
        ad_idx = find(getTDidx(trial_data,'epoch','ad'));
        for j = 1:length(pred_models)
            load(fullfile(rootDir,TDDir,['trainad_' pred_models{j} '_nokin'],['FF-PMd-M1_' filename '.mat']),'params');
            params.pred_array = pred_array;
            params.cov_array = cov_array;
            
            for unit = 1:size(neuron_id,1)
                neuron_idx = ismember(trial_data(1).([pred_array '_unit_guide']),neuron_id(unit,:),'rows');
                [y,x_full,~] = glm_prep_inputs(trial_data,unit,params.train_trials,params.train_start_idx,params.train_end_idx,params);
                b = glmfit(x_full,y,'poisson');
                
                for trial = ad_idx
                    [y,x_full,x_basic] = glm_prep_inputs(trial_data,unit,trial,params);
                    y = y(1:end-1); x_full = x_full(1:end-1,:);
                    yfit = exp([ones(size(x_full,1),1), x_full]*b)';
                    %                 bs = randi(length(yfit),length(yfit),1000);
                    %                 pr2 = prctile(compute_pseudo_R2(y(bs),yfit(bs),mean(y)),[2.5 97.5]);
                    %         rpr2 = prctile(compute_rel_pseudo_R2(y(bs),yfit_basic(bs),yfit(bs)),[2.5 97.5]);
                    trial_data(trial).Pred_spikes(:,unit) = yfit;
                    %                 trial_data(trial).pred_spikes_pr2(unit,:) = pr2;
                    %         trial_data(trial).pred_spikes_rpr2(unit,:) = rpr2;
                end
            end
            for trial = ad_idx
                trial_data(trial).(['Pred_' pred_models{j} '_spikes']) = trial_data(trial).Pred_spikes;
            end
        end
    end
    
    
    %%%%%%%%%%%%
    % PREDICT SPIKES
    load(fullfile(rootDir,TDDir,'trainad_nokin',['FF-PMd-M1_' filename '.mat']),'params');
    params.pred_array = pred_array;
    params.cov_array = cov_array;
    ad_idx = find(getTDidx(trial_data,'epoch','ad'));
    for unit = 1:size(neuron_id,1)
        tic;
        neuron_idx = ismember(trial_data(1).([pred_array '_unit_guide']),neuron_id(unit,:),'rows');
        [y,x_full,~] = glm_prep_inputs(trial_data,unit,params.train_trials,params.train_start_idx,params.train_end_idx,params);
        b = glmfit(x_full,y,'poisson');
        
        for trial = ad_idx
            [y,x_full,x_basic] = glm_prep_inputs(trial_data,unit,trial,params);
            y = y(1:end-1); x_full = x_full(1:end-1,:);
            yfit = exp([ones(size(x_full,1),1), x_full]*b)';
            %                 bs = randi(length(yfit),length(yfit),1000);
            %                 pr2 = prctile(compute_pseudo_R2(y(bs),yfit(bs),mean(y)),[2.5 97.5]);
            %         rpr2 = prctile(compute_rel_pseudo_R2(y(bs),yfit_basic(bs),yfit(bs)),[2.5 97.5]);
            trial_data(trial).Pred_spikes(:,unit) = yfit;
            %                 trial_data(trial).pred_spikes_pr2(unit,:) = pr2;
            %         trial_data(trial).pred_spikes_rpr2(unit,:) = rpr2;
        end
        toc;
    end
    for trial = ad_idx
        trial_data(trial).Pred_spikes = trial_data(trial).Pred_spikes;
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now build MIMO model from spikes to velocity in train trials
    input_vars = {'M1','Pred'};%'PMd','PCA','Potent','Null','Pred_potent','Pred_null'};
    output_vars = {'vel'};
    
    % get PCA/Potent/Null spaces
    %[V_potent, V_null, w_pmd, w_m1] = getPotentSpace(trial_data,'PMd','M1',params);
    
    % put together list of things to smooth with RCB?
    params.rcb_hpeaks = 0.05:0.05:0.5;
    params.rcb_b = 0.7;
    params.unit_lags = length(params.rcb_hpeaks)-1;
    
    ad_idx = find(getTDidx(trial_data,'epoch','ad'));
    ad_idx = ad_idx(1:floor(0.5*length(ad_idx)));
    
    temp_td = trial_data(params.train_trials);
    
%     for trial = 1:length(temp_td)
%         temp = sqrt(temp_td(trial).PMd_spikes);
%         temp = smoothSpikesForPCA(temp,params.bin_size,2*params.bin_size);
%         temp = temp*w_pmd;
%         temp_td(trial).PCA_spikes = temp(:,params.pca_dims.PMd);
%         temp_td(trial).Potent_spikes = temp(:,params.pca_dims.PMd)*V_potent;
%         temp_td(trial).Null_spikes = temp(:,params.pca_dims.PMd)*V_null;
%     end
%     
%     temp_td = convBasisFunc(temp_td,{'M1_spikes','PMd_spikes','PCA_spikes','Potent_spikes','Null_spikes'},params);
    
    temp_td = convBasisFunc(temp_td,{'M1_spikes','Pred_spikes'},params);

    for trial = 1:length(temp_td)
        t_idx = temp_td(trial).idx_go_cue - 2 : temp_td(trial).idx_trial_end -3;
        temp_td(trial).vel = temp_td(trial).vel(t_idx, :);
        temp_td(trial).M1_spikes = temp_td(trial).M1_spikes(t_idx, :);
        temp_td(trial).Pred_spikes = temp_td(trial).Pred_spikes(t_idx, :);
%         temp_td(trial).PMd_spikes = temp_td(trial).PMd_spikes(t_idx, :);
%         temp_td(trial).PCA_spikes = temp_td(trial).PCA_spikes(t_idx, :);
%         temp_td(trial).Potent_spikes = temp_td(trial).Potent_spikes(t_idx, :);
%         temp_td(trial).Null_spikes = temp_td(trial).Null_spikes(t_idx, :);
temp_td(trial).M1_spikes_shift = temp_td(trial).M1_spikes_shift(t_idx, :);
        temp_td(trial).Pred_spikes_shift = temp_td(trial).Pred_spikes_shift(t_idx, :);
%         temp_td(trial).PMd_spikes_shift = temp_td(trial).PMd_spikes_shift(t_idx, :);
%         temp_td(trial).PCA_spikes_shift = temp_td(trial).PCA_spikes_shift(t_idx, :);
%         temp_td(trial).Potent_spikes_shift = temp_td(trial).Potent_spikes_shift(t_idx, :);
%         temp_td(trial).Null_spikes_shift = temp_td(trial).Null_spikes_shift(t_idx, :);
    end
    y = cat(1,temp_td.vel);
    
    for j = 1:length(input_vars)
        x = [cat(1,temp_td.([input_vars{j} '_spikes'])), cat(1,temp_td.([input_vars{j} '_spikes_shift']))];
        
        % build wiener filter relating neurons to velocity
        [H,v,mcc]=filMIMO4(x,y,1,1,1);
        results(idx_file).(input_vars{j}).H = H;
        
        % cross-validate
        cv_folds = 10;
        folds = reshape(randperm(size(y,1)-rem(size(y,1),cv_folds)),[],cv_folds)';
        [cv_vaf, cv_r2] = deal(zeros(cv_folds,size(y,2)));
        for fold = 1:cv_folds
            test_idx = folds(fold,:);
            train_idx = true(1,size(y,1)); train_idx(test_idx) = false;
            y_train = y(train_idx,:); x_train = x(train_idx,:);
            x_test = x(test_idx,:); y_test = y(test_idx,:);
            [H,~,~]=filMIMO4(x_train,y_train,1,1,1);
            [pred,~,out_new] = predMIMO4(x_test,H,1,1,y_test);
            cv_vaf(fold,:) = calc_vaf(pred,out_new);
            for i = 1:size(pred,2)
                cv_r2(fold,i) = CalculateR2(pred(:,i),out_new(:,i));
            end
        end
        results(idx_file).(input_vars{j}).cv_vaf = cv_vaf;
        results(idx_file).(input_vars{j}).cv_r2 = cv_r2;
    end
    
    %%
    % Now get testing stuff
    temp_td = trial_data(ad_idx);
%     for trial = 1:length(temp_td)
%         temp = sqrt(temp_td(trial).PMd_spikes);
%         temp = smoothSpikesForPCA(temp,params.bin_size,2*params.bin_size);
%         temp = temp*w_pmd;
%         temp_td(trial).PCA_spikes = temp(:,params.pca_dims.PMd);
%         temp_td(trial).Potent_spikes = temp(:,params.pca_dims.PMd)*V_potent;
%         temp_td(trial).Null_spikes = temp(:,params.pca_dims.PMd)*V_null;
%     end

    temp_td = convBasisFunc(temp_td,{'M1_spikes','Pred_spikes'},params);
    for trial = 1:length(temp_td)
        t_idx = temp_td(trial).idx_go_cue - 2 : temp_td(trial).idx_trial_end -3;
        temp_td(trial).vel = temp_td(trial).vel(t_idx, :);
        temp_td(trial).M1_spikes = temp_td(trial).M1_spikes(t_idx, :);
        temp_td(trial).Pred_spikes = temp_td(trial).Pred_spikes(t_idx, :);
%         temp_td(trial).PMd_spikes = temp_td(trial).PMd_spikes(t_idx, :);
%         temp_td(trial).PCA_spikes = temp_td(trial).PCA_spikes(t_idx, :);
%         temp_td(trial).Potent_spikes = temp_td(trial).Potent_spikes(t_idx, :);
%         temp_td(trial).Null_spikes = temp_td(trial).Null_spikes(t_idx, :);
        temp_td(trial).M1_spikes_shift = temp_td(trial).M1_spikes_shift(t_idx, :);
        temp_td(trial).Pred_spikes_shift = temp_td(trial).Pred_spikes_shift(t_idx, :);
%         temp_td(trial).PMd_spikes_shift = temp_td(trial).PMd_spikes_shift(t_idx, :);
%         temp_td(trial).PCA_spikes_shift = temp_td(trial).PCA_spikes_shift(t_idx, :);
%         temp_td(trial).Potent_spikes_shift = temp_td(trial).Potent_spikes_shift(t_idx, :);
%         temp_td(trial).Null_spikes_shift = temp_td(trial).Null_spikes_shift(t_idx, :);
    end
    
    bins = 1:5:length(temp_td);
    
    for j = 1:length(input_vars)
        vaf = zeros(length(bins)-1,size(y,2));
        r2 = zeros(length(bins)-1,size(y,2));
        for bin = 1:length(bins)-1
            trial = bins(bin):bins(bin+1)-1;
            y = cat(1,temp_td(trial).vel);
            
            x = [cat(1,temp_td(trial).([input_vars{j} '_spikes'])), cat(1,temp_td(trial).([input_vars{j} '_spikes_shift']))];
            
            [pred,in_new,out_new] = predMIMO4(x,results(idx_file).(input_vars{j}).H,1,1,y);
            
            % get r2
            vaf(bin,:) = calc_vaf(pred,out_new);
            for i = 1:size(pred,2)
                r2(bin,i) = CalculateR2(pred(:,i),out_new(:,i));
            end
        end
        results(idx_file).(input_vars{j}).r2 = r2;
        results(idx_file).(input_vars{j}).vaf = vaf;
    end
    
end

%% Plot
figure;
subplot1(1,length(input_vars));
for j = 1:length(input_vars)
    input_var = input_vars{j};
    temp = zeros(1,length(results));
    for i = 1:length(results)
        temp(i) = size(results(i).(input_var).r2,1);
    end
    
    n_bins = min(temp);
    d = zeros(n_bins,size(results(1).(input_var).r2,2),length(results));
    for i = 1:length(results)
        d(:,:,i) = (results(i).(input_var).vaf(1:n_bins,:));% - repmat(mean(results(i).(input_var).cv_r2,1),n_bins,1))./repmat(mean(results(i).(input_var).cv_r2,1),n_bins,1);
    end
    
    h(j) = subplot1(j);
    hold all;
    m = mean(d(:,1,:),3);
    s = std(d(:,1,:),[],3)./sqrt(size(d,3));
    plot(m,'r-','LineWidth',3);
    plot(m-s,'r--','LineWidth',2);
    plot(m+s,'r--','LineWidth',2);
    
    m = mean(d(:,2,:),3);
    s = std(d(:,2,:),[],3)./sqrt(size(d,3));
    plot(m,'b-','LineWidth',3);
    plot(m-s,'b--','LineWidth',2);
    plot(m+s,'b--','LineWidth',2);
    
    set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[0 1]);
    title(input_var);
end
linkaxes(h);






