% This will plot mean prediction R2 for the first and last trials of some
% epoch for different models for all neurons pooled
dataSummary;

sessions = { ...
    'Chewie','2016-09-15'; ... % CF
    'Chewie','2016-10-05'; ...
    'Chewie','2016-10-07'; ...
    'Chewie','2016-10-11'; ...
%     'Mihili','2014-02-03'; ...
%     'Mihili','2014-02-17'; ...
%     'Mihili','2014-02-18'; ...
%     'Mihili','2014-03-07'; ...
    'Chewie','2016-09-09'; ... % VR
    'Chewie','2016-09-12'; ...
    'Chewie','2016-09-14'; ...
    'Chewie','2016-10-06'; ...
    'Mihili','2014-03-03'; ...
    'Mihili','2014-03-04'; ...
    'Mihili','2014-03-06'; ...
    };

basenames = {'trainad','trainad'};
extranames = {'singleneurons','singleneurons'};
array_models = {'M1-M1','PMd-M1'};


pert = 'FF';
tasks = {'CO'};
dates = sessions(:,2);
monkeys = unique(sessions(:,1));


which_metric = 'rpr2'; % 'rpr2','pr2_full','pr2_basic'
pr2_cutoff = 0.0;
pr2_op = 'min'; % which operation for filtering ('min','max','mean','median')
do_same_cells = false;
basic_pr2_check = false;


test_trials = {'AD',[0 0.5]};
num_trials = 20;


%%
plot_colors = [0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    %     0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840];

session_idx = ismember(filedb.Monkey,monkeys) & ismember(filedb.Perturbation,pert) & ismember(filedb.Task,tasks) & ismember(filedb.Date,dates);


% if we want the same cells for all models, make sure all are good
if do_same_cells
    all_good_cells = cell(1,length(basenames));
    for idx_cond = 1:length(basenames)
        if isempty(extranames{idx_cond})
            outputSubdir = basenames{idx_cond};
        else
            outputSubdir = [basenames{idx_cond} '_' extranames{idx_cond}];
        end
        
        idx = find(strcmpi(pert,filedb.Perturbation) & session_idx);
        filenames = cell(1,length(idx));
        for s = 1:length(idx)
            filenames{s} = [filedb.Monkey{idx(s)} '_' filedb.Task{idx(s)} '_' filedb.Perturbation{idx(s)} '_' filedb.Date{idx(s)}];
        end
        
        % build list of filenames
        good_cells = [];
        for file = 1:length(filenames)
            % get good cells for both conditions
            out_struct = get_plot_metrics({fullfile(rootDir,resultsDir,outputSubdir,[pert '-' array_models{idx_cond} '_' filenames{file} '.mat'])}, ...
                struct( ...
                'which_metric',which_metric, ...
                'epochs',{epochs}, ...
                'pr2_cutoff',pr2_cutoff, ...
                'pr2_op',pr2_op, ...
                'pr2_ad_check', pr2_ad_check, ...
                'do_good_cells',true, ...
                'do_behavior',false, ...
                'filter_trials',filter_trials, ...
                'basic_pr2_check',basic_pr2_check));
            
            good_cells = [good_cells; out_struct.good_cells];
        end
        all_good_cells{idx_cond} = good_cells;
    end
    good_cells = all(cell2mat(all_good_cells),2);
end


for idx_cond = 1:length(basenames)
    if isempty(extranames{idx_cond})
        outputSubdir = basenames{idx_cond};
    else
        outputSubdir = [basenames{idx_cond} '_' extranames{idx_cond}];
    end
    
    idx = find(strcmpi(pert,filedb.Perturbation) & session_idx);
    filenames = cell(1,length(idx));
    for s = 1:length(idx)
        filenames{s} = [filedb.Monkey{idx(s)} '_' filedb.Task{idx(s)} '_' filedb.Perturbation{idx(s)} '_' filedb.Date{idx(s)}];
    end
    
    % for each good cell, predict in the first NUM_TRIALS of EPOCH and last
    % NUM_TRIALS of EPOCH
    all_r2 = [];
    for file = 1:length(filenames)
        filepath = fullfile(rootDir,resultsDir,outputSubdir,[pert '-' array_models{idx_cond} '_' filenames{file} '.mat']);
        out_struct = get_plot_metrics({filepath}, ...
            struct( ...
            'which_metric',which_metric, ...
            'epochs',{test_trials(1)}, ...
            'pr2_cutoff',pr2_cutoff, ...
            'pr2_op',pr2_op, ...
            'pr2_ad_check', false, ...
            'do_good_cells',~do_same_cells, ...
            'do_behavior',false, ...
            'filter_trials',false, ...
                'basic_pr2_check',basic_pr2_check));
        
        params = out_struct.params;
        
        if do_same_cells
            cv = out_struct.cv(good_cells,:);
            total_significant = sum(good_cells);
            total_cells = length(good_cells);
        else
            cv = out_struct.cv;
            total_significant = out_struct.total_significant;
            total_cells = out_struct.total_cells;
        end
        
        load(fullfile(rootDir,resultsDir,outputSubdir,[pert '-' array_models{idx_cond} '_' filenames{file} '_td.mat']));
        
        % e.g. trained on 0.5->1 of AD
        trial_idx = getTDidx(trial_data_test,'epoch',test_trials{1},'range',test_trials{2});
        
        y = cat(1,trial_data_test(trial_idx(1:num_trials)).y);
        yb = cat(1,trial_data_test(trial_idx(1:num_trials)).yfit_basic);
        yf = cat(1,trial_data_test(trial_idx(1:num_trials)).yfit_full);
        
        if do_same_cells
            y = y(:,good_cells);
            yb = yb(:,good_cells);
            yf = yf(:,good_cells);
        end
        
        r2 = zeros(size(y,2),2);
        for unit = 1:size(y,2)
            r2(unit,1) = compute_rel_pseudo_R2(y(:,unit),yb(:,unit),yf(:,unit));
        end
        
        y = cat(1,trial_data_test(trial_idx(end-num_trials:end)).y);
        yb = cat(1,trial_data_test(trial_idx(end-num_trials:end)).yfit_basic);
        yf = cat(1,trial_data_test(trial_idx(end-num_trials:end)).yfit_full);
        
        if do_same_cells
            y = y(:,good_cells);
            yb = yb(:,good_cells);
            yf = yf(:,good_cells);
        end
        
        for unit = 1:size(y,2)
            r2(unit,2) = compute_rel_pseudo_R2(y(:,unit),yb(:,unit),yf(:,unit));
        end
        
        all_r2 = [all_r2; r2];
    end
    
    figure;
    subplot(2,1,1);
    hist(all_r2(:,1),-1:0.1:0.5);
    axis('tight');
    set(gca,'Box','off','TickDir','out','FontSize',14);
    
    subplot(2,1,2);
    hist(all_r2(:,2),-1:0.1:0.5);
    axis('tight');
    set(gca,'Box','off','TickDir','out','FontSize',14);
    
    xlabel(which_metric);
    title(outputSubdir);
end








