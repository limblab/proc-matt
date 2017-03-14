% This will plot example predictions for single neurons for the full/basic
% models
dataSummary;

sessions = { ...
    'Chewie','2016-10-07'; ...
    };

basename = 'trainad';
extraname = 'potent_bl';
array_model = 'PMd-M1';

pert = 'FF';
tasks = {'CO'};
dates = sessions(:,2);
monkeys = unique(sessions(:,1));


which_metric = 'rpr2'; % 'rpr2','pr2_full','pr2_basic'
pr2_cutoff = 0.01;
pr2_op = 'min'; % which operation for filtering ('min','max','mean','median')
do_same_cells = false;

epoch = {'AD'};
num_trials = 20;

%%
font_size = 14;
line_width = 2;

plot_colors = [0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    %     0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840];

session_idx = ismember(filedb.Monkey,monkeys) & ismember(filedb.Perturbation,pert) & ismember(filedb.Task,tasks) & ismember(filedb.Date,dates);

% if we want the same cells for all models, make sure all are good
same_good_cells;

if isempty(extraname)
    outputSubdir = basename;
else
    outputSubdir = [basename '_' extraname];
end

idx = find(strcmpi(pert,filedb.Perturbation) & session_idx);

filename = [filedb.Monkey{idx} '_' filedb.Task{idx} '_' filedb.Perturbation{idx} '_' filedb.Date{idx}];

% for each good cell, predict in the first NUM_TRIALS of EPOCH and last
% NUM_TRIALS of EPOCH
all_r2 = [];

filepath = fullfile(rootDir,resultsDir,outputSubdir,[pert '-' array_model '_' filename '.mat']);
out_struct = get_plot_metrics(filepath, ...
    struct( ...
    'which_metric',which_metric, ...
    'epochs',{epoch}, ...
    'pr2_cutoff',pr2_cutoff, ...
    'pr2_op',pr2_op, ...
    'pr2_ad_check', false, ...
    'do_good_cells',~do_same_cells, ...
    'do_behavior',false, ...
    'filter_trials',false));

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

load(fullfile(rootDir,resultsDir,outputSubdir,[pert '-' array_model '_' filename '_td.mat']));
% trained on 0.5->1 of AD
ad_idx = getTDidx(trial_data_test,'epoch','AD',[0 0.5]);

for unit = 1:size(cv,2)
    figure;
    r2 = zeros(3,2);
    
    y1 = cat(1,trial_data_test(ad_idx(1:num_trials)).([pred_array '_spikes']));
    yf1 = cat(1,trial_data_test(ad_idx(1:num_trials)).yfit_full);
    yb1 = cat(1,trial_data_test(ad_idx(1:num_trials)).yfit_basic);
    
    if do_same_cells
        y1 = y1(:,good_cells);
        yf1 = yf1(:,good_cells);
        yb1 = yb1(:,good_cells);
    end
    
    y1 = y1(:,unit);
    yb1 = yb1(:,unit);
    yf1 = yf1(:,unit);
    
    r2(1,1) = compute_pseudo_r2(y1,yb1,mean(y1));
    r2(2,1) = compute_pseudo_r2(y1,yf1,mean(y1));
    r2(3,1) = compute_rel_pseudo_r2(y1,yb1,yf1);
    
    y2 = cat(1,trial_data_test(ad_idx(end-num_trials:end)).([pred_array '_spikes']));
    yf2 = cat(1,trial_data_test(ad_idx(end-num_trials:end)).yfit_full);
    yb2 = cat(1,trial_data_test(ad_idx(end-num_trials:end)).yfit_basic);
    
    if do_same_cells
        y2 = y2(:,good_cells);
        yf2 = yf2(:,good_cells);
        yb2 = yb2(:,good_cells);
    end
    y2 = y2(:,unit);
    yb2 = yb2(:,unit);
    yf2 = yf2(:,unit);
    
    r2(1,2) = compute_pseudo_r2(y2,yb2,mean(y2));
    r2(2,2) = compute_pseudo_r2(y2,yf2,mean(y2));
    r2(3,2) = compute_rel_pseudo_r2(y2,yb2,yf2);
    
    max_y = max([y1 yb1 yf1 y2 yb2 yf2]);
    
    subplot(1,2,1); hold all;
    plot(y1,'k','LineWidth',line_width);
    plot(yb1,'b','LineWidth',line_width);
    plot(yf1,'r','LineWidth',line_width);
    set(gca,'Box','off','TickDir','out','FontSize',font_size,'YLim',[0 max_y]);
    text(num_trials/2,max_y-1,['pR2b = ' num2str(r2(1,1))]);
    text(num_trials/2,max_y-1,['pR2f = ' num2str(r2(2,1))]);
    text(num_trials/2,max_y-1,['rpR2 = ' num2str(r2(3,1))]);
    
    
    subplot(1,2,2); hold all;
    plot(y2,'k','LineWidth',line_width);
    plot(yb2,'b','LineWidth',line_width);
    plot(yf2,'r','LineWidth',line_width);
    set(gca,'Box','off','TickDir','out','FontSize',font_size,'YLim',[0 max_y]);
    text(num_trials/2,max_y-1,['pR2b = ' num2str(r2(1,2))]);
    text(num_trials/2,max_y-1,['pR2f = ' num2str(r2(2,2))]);
    text(num_trials/2,max_y-1,['rpR2 = ' num2str(r2(3,2))]);
    
    pause;
    close all;
    
end





