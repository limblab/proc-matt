%%
clear
clc;
close all;

if 0
    a = reshape(e_pr2{2},numel(e_pr2{2}),1);
    g = reshape( repmat(1:size(e_pr2{2},2),size(e_pr2{2},1),1), numel(e_pr2{2}), 1 );
    anovan(a,g,'Display','off')
end

dataSummary;

% % M1PMd sessions
sessions = { ...
    'Chewie','2016-09-09'; ... % VR
    'Chewie','2016-09-12'; ...
    'Chewie','2016-09-14'; ...
    'Chewie','2016-10-06'; ...
    'Mihili','2014-03-03'; ...
    'Mihili','2014-03-04'; ...
    'Mihili','2014-03-06'; ...
    'Mihili','2015-06-23'; ...
    'Mihili','2015-06-25'; ...
    'Mihili','2015-06-26'; ...
    'Chewie','2016-09-15'; ... % CF
    'Chewie','2016-10-05'; ...
    'Chewie','2016-10-07'; ...
    'Chewie','2016-10-11'; ...
    'Mihili','2014-02-03'; ...
    'Mihili','2014-02-17'; ...
    'Mihili','2014-02-18'; ...
    'Mihili','2014-03-07'; ...
    'Mihili','2015-06-15'; ...
    'Mihili','2015-06-16'; ...
    'Mihili','2015-06-17'; ...
        'MrT','2013-08-19'; ... % CF
        'MrT','2013-08-21'; ...
        'MrT','2013-08-23'; ...
        'MrT','2013-09-03'; ... %VR
        'MrT','2013-09-05'; ...
        'MrT','2013-09-09'; ...
%     };
% % CHEWIE M1 ONLY
% % sessions = { ...
    'Chewie','2013-10-22'; ...
    'Chewie','2013-10-23'; ...
    'Chewie','2013-10-31'; ...
    'Chewie','2013-11-01'; ...
    'Chewie','2013-12-03'; ...
    'Chewie','2013-12-20'; ...
         'Chewie','2015-06-29'; ... % short washout
    'Chewie','2015-06-30'; ...
    'Chewie','2015-07-01'; ...
    'Chewie','2015-07-03'; ...
    'Chewie','2015-07-06'; ...
          'Chewie','2015-07-13'; ... % short washout
    'Chewie','2013-10-03'; ...
    'Chewie','2015-07-10'; ...
    'Chewie','2015-07-14'; ...
    'Chewie','2015-07-15'; ...
    'Chewie','2015-07-16'; ...
    };




basenames = {'trainad','trainad','trainad','trainad','trainad'};
extranames = {'','','','potent_bl','null_bl'};
array_models = {'M1-M1','PMd-PMd','PMd-M1','PMd-M1','PMd-M1'};

% the good potent/null sessions
sessions = { ...
    'Chewie','2016-09-15'; ... % CF
    'Chewie','2016-10-05'; ...
    'Chewie','2016-10-07'; ...
    'Chewie','2016-10-11'; ...
    'Mihili','2014-02-03'; ...
    'Mihili','2014-02-17'; ...
    'Mihili','2014-02-18'; ...
    'Mihili','2014-03-07'; ...
%         'Chewie','2016-09-09'; ... % VR
%         'Chewie','2016-09-12'; ...
%         'Chewie','2016-09-14'; ...
%         'Chewie','2016-10-06'; ...
%         'Mihili','2014-03-03'; ...
%         'Mihili','2014-03-04'; ...
%         'Mihili','2014-03-06'; ...
    };
% basenames = {'trainad','trainad','trainad'};
% extranames = {'','potent_bl','null_bl'};
% array_models = {'PMd-M1','PMd-M1','PMd-M1'};
% basenames = {'trainad'};
% extranames = {'potent_null_pred'};
% array_models = {'null-potent'};

pert = 'FF';
tasks = {'CO'};
dates = sessions(:,2);
monkeys = unique(sessions(:,1));

which_metric = 'pr2_basic'; % 'rpr2','pr2_full','pr2_basic'
pr2_cutoff = 0;
pr2_op = 'min'; % which operation for filtering ('min','max','mean','median')
pr2_ad_check = false; % only keeps cells that predict in WO
do_same_cells = false; % really only works for testing tweaks of model

plot_op = 'mean';
group_size = 20;
how_to_group = 'slide'; % average, pool, slide

do_norm = true;
do_diff = true;
remove_outliers = true;
error_bars = 'ste'; % 'boot','ste'
num_bootstraps = 1000;

do_regression_line = false;
add_plot = true;
do_subplot = false;
filter_trials = false;

epochs = {'AD'};%{'BL','AD','WO'};

if do_norm
    min_y = -0.7; max_y = 0.1;
else
    min_y = -0.1; max_y = 0.1;
end

%%

plot_colors = [0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    %     0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840];

session_idx = ismember(filedb.Monkey,monkeys) & ismember(filedb.Perturbation,pert) & ismember(filedb.Task,tasks) & ismember(filedb.Date,dates);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                'filter_trials',filter_trials));
            
            good_cells = [good_cells; out_struct.good_cells];
        end
        all_good_cells{idx_cond} = good_cells;
    end
    good_cells = all(cell2mat(all_good_cells),2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold all;
h = zeros(1,length(basenames)*2);
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
    filepaths = cell(1,length(filenames));
    for file = 1:length(filenames)
        filepaths{file} = fullfile(rootDir,resultsDir,outputSubdir,[pert '-' array_models{idx_cond} '_' filenames{file} '.mat']);
    end
    
    out_struct = get_plot_metrics(filepaths, ...
        struct( ...
        'which_metric',which_metric, ...
        'epochs',{epochs}, ...
        'pr2_cutoff',pr2_cutoff, ...
        'pr2_op',pr2_op, ...
        'pr2_ad_check', pr2_ad_check, ...
        'do_good_cells',~do_same_cells, ...
        'do_behavior',false, ...
        'filter_trials',filter_trials));
    
    if do_same_cells
        cv = out_struct.cv(good_cells,:);
        e_pr2 = cell(1,length(out_struct.e_pr2));
        for i = 1:length(out_struct.e_pr2)
            e_pr2{i} = out_struct.e_pr2{i}(good_cells,:);
        end
        e_inds = out_struct.e_inds;
        total_significant = sum(good_cells);
        total_cells = length(good_cells);
    else
        cv = out_struct.cv;
        e_pr2 = out_struct.e_pr2;
        e_inds = out_struct.e_inds;
        total_significant = out_struct.total_significant;
        total_cells = out_struct.total_cells;
    end
    
    disp([outputSubdir ' - ' array_models{idx_cond} ' - % cells with significant rel-pseudo-r2: ' num2str(total_significant) '/' num2str(total_cells)]);
    disp([outputSubdir ' - ' array_models{idx_cond} ' - Mean baseline metric: ' num2str(mean(nanmean(cv,2))) ' +/- ' num2str(std(nanmean(cv,2)))]);
    
    subplot(1,length(basenames),idx_cond); hold all;
    hist(mean(cv,2),-1:0.1:1);
    axis('tight');
    V = axis;
    plot([0.025 0.025],V(3:4),'r-','LineWidth',1);
    set(gca,'Box','off','TickDir','out','FontSize',14);
    xlabel(which_metric);
    ylabel('Count');
    title(outputSubdir);
end