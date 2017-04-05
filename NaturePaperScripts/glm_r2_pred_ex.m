close all; clear; clc;
dataSummary;

% the good potent/null sessions
sessions = { ...
    %     'Chewie','2016-09-15'; ... % CF
    %     'Chewie','2016-09-21'; ...
    %     'Chewie','2016-10-05'; ...
    'Chewie','2016-10-07'; ...
    %     'Chewie','2016-10-11'; ...
    %     'Mihili','2014-02-03'; ... 
    %     'Mihili','2014-02-17'; ...
    %     'Mihili','2014-02-18'; ...
    %     'Mihili','2014-03-07'; ...
    %     'Chewie','2016-09-09'; ... % VR
    %     'Chewie','2016-09-12'; ...
    %     'Chewie','2016-09-14'; ...
    %     'Chewie','2016-10-06'; ...
    %     'Mihili','2014-03-03'; ...
    %     'Mihili','2014-03-04'; ...
    %     'Mihili','2014-03-06'; ...
    };

basenames = {'trainad'};
extranames = {'null8'};
array_models = {'PMd-M1'};

pert            = 'FF';
tasks           = {'CO'};
dates           = sessions(:,2);
monkeys         = unique(sessions(:,1));

% last one is "cross validated" one
test_trials = {'AD',[1,10], [71,80]};
num_bootstraps = 0;

which_metric    = 'rpr2'; % 'rpr2','pr2_full','pr2_basic'
min_r2      =  0.05;

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
        
        
        % loop along neurons
        for unit = 1:length(results.bl_model)
            % set up the figure, subplot for each block
            figure;
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
                
                temp1 = y(:,unit);
                temp2 = yh_b(:,unit);
                temp3 = yh_f(:,unit);
                
                pr2_b = mean(compute_pseudo_R2(temp1(bs),temp2(bs),mean(temp1)));
                pr2_f = mean(compute_pseudo_R2(temp1(bs),temp3(bs),mean(temp1)));
                rpr2 = mean(compute_rel_pseudo_R2(temp1(bs),temp2(bs),temp3(bs)));
                
                % do the plotting
                subplot(length(test_trials)-1,1,i); hold all;
                plot(temp1,'k','LineWidth',2)
                plot(temp2,'LineWidth',2);
                plot(temp3,'LineWidth',2);
                set(gca,'Box','off','TickDir','out','FontSize',14);
                legend( {...
                    'Actual', ...
                    ['Basic (pR2=' num2str(pr2_b) ')'], ...
                    ['Full (pr2=' num2str(pr2_f) '; rpR2=' num2str(rpr2) ')']});
                if i == 1
                    temp = trial_data_test(1).([params.pred_array '_unit_guide']);
                    title(temp(unit,:)')
                end
            end
            if rpr2 > min_r2
                pause;
            end
            close all;
        end
    end
end