close all; clear; clc;
dataSummary;

% the good potent/null sessions
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

basenames = {'trainad'};
extranames = {'singleneurons'};
array_models = {'M1-M1'};

do_the_plots = false; % will pot neuron examples

pert            = 'FF';
tasks           = {'CO'};
dates           = sessions(:,2);
monkeys         = unique(sessions(:,1));

% last one is "cross validated" one
test_trials = {'AD',[1,10], [71,80]};
num_bootstraps = 0;

which_metric    = 'rpr2'; % 'rpr2','pr2_full','pr2_basic'
pr2_cutoff      =  0.0;
pr2_op          = 'min'; % which operation for filtering ('min','max','mean','median')
min_fr          = 0.15;
basic_pr2_check = false;

n_bins = 15;
do_diff = false;
do_norm = false;

rem_outliers = false;
num_std = 10;
only_high_fr_cells = false;
only_good_late_cells = false;
plot_diff = false;

session_idx = ismember(filedb.Monkey,monkeys) & ismember(filedb.Perturbation,pert) & ismember(filedb.Task,tasks) & ismember(filedb.Date,dates);


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
    filepaths = cell(1,length(filenames));
    for file = 1:length(filenames)
        filepaths{file} = fullfile(rootDir,resultsDir,outputSubdir,[pert '-' array_models{idx_cond} '_' filenames{file} '.mat']);
    end
    
    [file_r2,file_fr,file_r2_basic] = deal(cell(length(filenames),1));
    for file = 1:length(filenames)
        load(fullfile(rootDir,resultsDir,outputSubdir,[pert '-' array_models{idx_cond} '_' filenames{file} '.mat']),'results','params');
        load(fullfile(rootDir,resultsDir,outputSubdir,[pert '-' array_models{idx_cond} '_' filenames{file} '_td.mat']),'trial_data_test');
        
        [r2,fr,r2_basic] = deal(zeros(length(results.bl_model),length(test_trials)-1));
        for i = 1:length(test_trials)-1
            [~,td] = getTDidx(trial_data_test,'epoch',test_trials{1},'range',test_trials{i+1});
            
            y = cat(1,td.y);
            yh_b = cat(1,td.yfit_basic);
            switch lower(which_metric)
                case 'pr2_full'
                    yh_f = cat(1,td.yfit_full);
                case 'rpr2'
                    yh_f = cat(1,td.yfit_full);
            end
            
            if num_bootstraps < 2
                bs = 1:size(y,1);
            else
                bs = randi(size(y,1),size(y,1),num_bootstraps);
            end
            
            for unit = 1:length(results.bl_model)
                
                fr(unit,i) = mean(y(:,unit));
                
                temp1 = y(:,unit);
                temp2 = yh_b(:,unit);
                temp3 = yh_f(:,unit);
                r2_basic(unit,i) = mean(compute_pseudo_R2(temp1(bs),temp2(bs),mean(temp1)));
                switch lower(which_metric)
                    case 'pr2_basic'
                        r2(unit,i) = mean(compute_pseudo_R2(temp1(bs),temp2(bs),mean(temp1)));
                    case 'pr2_full'
                        r2(unit,i) = mean(compute_pseudo_R2(temp1(bs),temp3(bs),mean(temp1)));
                    case 'rpr2'
                        r2(unit,i) = mean(compute_rel_pseudo_R2(temp1(bs),temp2(bs),temp3(bs)));
                end
                
                if do_the_plots
                    figure; hold all;
                    plot(temp1,'k','LineWidth',2)
                    plot(temp2,'LineWidth',2);
                    plot(temp3,'LineWidth',2);
                    set(gca,'Box','off','TickDir','out','FontSize',14);
                    legend( {...
                        'Actual', ...
                        ['Basic (pR2=' num2str(mean(compute_pseudo_R2(temp1(bs),temp2(bs),mean(temp1)))) ')'], ...
                        ['Full (rpR2=' num2str(mean(compute_rel_pseudo_R2(temp1(bs),temp2(bs),temp3(bs)))) ')']});
                    pause;
                    close all;
                end
            end
        end
        file_r2{file} = r2;
        file_fr{file} = fr;
        file_r2_basic{file} = r2_basic;
    end
    
    r2_b = cell2mat(file_r2_basic);
    fr = cell2mat(file_fr);
    r2 = cell2mat(file_r2);
    r2_e = r2(:,1);
    r2_l = r2(:,2);
    
    out_struct = get_plot_metrics(filepaths, ...
        struct( ...
        'which_metric',which_metric, ...
        'epochs',{test_trials(1)}, ...
        'pr2_cutoff',pr2_cutoff, ...
        'pr2_op',pr2_op, ...
        'pr2_ad_check', false, ...
        'do_good_cells',false, ...
        'do_behavior',false, ...
        'filter_trials',false, ...
        'basic_pr2_check',basic_pr2_check));
    good_cells = out_struct.good_cells > 0;
    
    if only_high_fr_cells
        good_cells = good_cells & all(fr > min_fr,2);
    end
    
    if only_good_late_cells
        good_cells = good_cells & r2(:,end) > 0;
    end
    
    good_cells = good_cells & all(r2_b > 0,2);
    
    cv = mean(out_struct.cv,2);
    
    if do_diff
        r2_e = r2_e - repmat(cv,1,size(r2_e,2));
        r2_l = r2_l - repmat(cv,1,size(r2_l,2));
    end
    if do_norm
        r2_e = r2_e ./ repmat(cv,1,size(r2_e,2));
        r2_l = r2_l ./ repmat(cv,1,size(r2_l,2));
    end
    
    r2_e(abs(r2_e) == Inf) = NaN;
    r2_l(abs(r2_l) == Inf) = NaN;
    
    if rem_outliers
        r2_e(abs(r2_e) > 20) = NaN;
        r2_l(abs(r2_l) > 20) = NaN;
        
        bad_idx = r2_e(:) > num_std*nanstd(r2_e(good_cells));
        bad_idx = bad_idx | r2_l > num_std*nanstd(r2_l(good_cells));
        good_cells = good_cells & ~bad_idx;
    end
    
    
    if plot_diff
        x_max = nanmax(r2_e(good_cells) - r2_l(good_cells));
        x_min = nanmin(r2_e(good_cells) - r2_l(good_cells));
    else
        x_max = max([nanmax(r2_e(good_cells)), nanmax(r2_l(good_cells))]);
        x_min = min([nanmin(r2_e(good_cells)), nanmin(r2_l(good_cells))]);
    end
    
    plot_out(idx_cond).r2_e = r2_e;
    plot_out(idx_cond).r2_l = r2_l;
    plot_out(idx_cond).cv = cv;
    plot_out(idx_cond).x_lim = [x_min,x_max];
    plot_out(idx_cond).good_cells = good_cells;
    
end

% plot now
x_min = Inf;
x_max = -Inf;
for i = 1:length(plot_out)
    x_min = min([x_min, plot_out(i).x_lim(1)]);
    x_max = max([x_max, plot_out(i).x_lim(2)]);
end
bins = x_min:(x_max-x_min)/n_bins:x_max;

%%
figure;
%subplot1(length(plot_out),1);
for i = 1:length(plot_out)
    r2_e = plot_out(i).r2_e;
    r2_l = plot_out(i).r2_l;
    
    %     r2_e(r2_e < -0.5) = NaN;
    %     r2_l(r2_l < -0.5) = NaN;
    %     x_min = -4; x_max = 8;
    %     bins = -4:0.4:8;
    cv = plot_out(i).cv;
    good_cells = plot_out(i).good_cells;
    
    %subplot1(i); hold all;
    subplot(length(plot_out),1,i); hold all;
    if ~plot_diff
        hist(r2_e(good_cells),bins);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','r');
        hist(r2_l(good_cells),bins);
        h = findobj(gca,'Type','patch');
        set(h,'EdgeColor','w','FaceAlpha',0.7,'EdgeAlpha',0.7);
        set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[x_min,x_max]);
        [~,p] = ttest2(r2_l(good_cells),r2_e(good_cells));
        title(p);
        
        if i == length(plot_out)
            legend({'First 10 Trials','Last 10 Trials'},'Location','NorthWest');
        end
    else
        hist((r2_e(good_cells) - r2_l(good_cells)),bins);
        h = findobj(gca,'Type','patch');
        set(h,'EdgeColor','w','FaceAlpha',0.7,'EdgeAlpha',0.7);
        set(gca,'Box','off','TickDir','out','FontSize',14);
        V = axis;
        plot([0,0],V(3:4),'k--');
        [~,p] = ttest((r2_e(good_cells)-r2_l(good_cells)),0,'tail','left');
        title(p);
    end
    ylabel(extranames{i});
    
    
end
xlabel('relative pseudo-R2');

