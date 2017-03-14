function out_struct = get_plot_metrics(filenames,inparams)

which_metric = inparams.which_metric;
epochs = inparams.epochs;
pr2_cutoff = inparams.pr2_cutoff;
pr2_op = inparams.pr2_op;
pr2_ad_check = inparams.pr2_ad_check;
do_good_cells = inparams.do_good_cells;
do_behavior = inparams.do_behavior;
filter_trials = inparams.filter_trials;

[cv,good_cells] = deal([]);
% loop along files and group together
[total_cells,total_significant] = deal(0);
[e_pr2,e_inds,behavior,good_trials] = deal(cell(1,length(epochs)));
for file = 1:length(filenames) % loop along sessions
    if exist(filenames{file},'file')
        load(filenames{file},'results','params');
        
        if ischar(params.test_epochs(1,:))
            a = cell(size(params.test_epochs,1),1);
            for i = 1:size(params.test_epochs,1), a{i} = params.test_epochs(i,:); end, params.test_epochs = a';
            params.test_epochs = a;
        end
        
        
        switch lower(pr2_op)
            case 'mean'
                temp_pr2 = mean(mean(results.([which_metric '_cv']),3),2);
            case 'median'
                temp_pr2 = median(mean(results.([which_metric '_cv']),3),2);
            case 'max'
                temp_pr2 = max(mean(results.([which_metric '_cv']),3),[],2);
            case 'min'
                temp_pr2 = min(mean(results.([which_metric '_cv']),3),[],2);
        end
        
        temp_metric = results.(which_metric);
        temp_metric_cv = results.([which_metric '_cv']);
        
        if any(any(any(temp_metric))) || any(any(any(temp_metric_cv)))
            disp('isnan...');
        end
        
        
        if pr2_ad_check % only take cells that can predict significantly in first bin of AD
            temp_idx = find(strcmpi(params.test_epochs,'AD'));
            good_idx = temp_metric(:,temp_idx(end),1) > pr2_cutoff;
        else % take the ones that are significant from cross validation
            good_idx = temp_pr2(:,1) > pr2_cutoff & min(mean(results.pr2_full_cv,3),[],2) > pr2_cutoff & min(mean(results.pr2_basic_cv,3),[],2) > pr2_cutoff ;
        end
        
        if do_good_cells
            temp_metric = temp_metric(good_idx,:,:);
            temp_metric_cv = temp_metric_cv(good_idx,:,:);
        end
        
        % get some stats
        total_significant = total_significant+sum(good_idx);
        total_cells = total_cells+length(good_idx);
        
        
        
        for e = 1:length(epochs)
            e_inds{e} = find(strcmpi(params.test_epochs,epochs{e}));
            if file == 1
                e_pr2{e} = mean(temp_metric(:,strcmpi(params.test_epochs,epochs{e}),:),3);
            else
                e_inds{e} = e_inds{e}(1:min([length(e_inds{e}), size(e_pr2{e},2)]));
                e_pr2{e} = [e_pr2{e}(:,1:length(e_inds{e})); mean(temp_metric(:,e_inds{e},:),3)];
            end
        end
        
        cv = [cv; mean(temp_metric_cv,3)];
        good_cells = [good_cells; good_idx];
        
        
        if do_behavior
            % this is super hacky but it works. Fix it after SfN
            [temp_path, temp_name, temp_ext]=fileparts(filenames{file});
            temp_path = strsplit(temp_path,filesep);
            new_path = temp_path{1};
            for i = 2:length(temp_path)-1
                new_path = [new_path, filesep, temp_path{i}];
            end
            load(fullfile(new_path,[temp_name(11:end) temp_ext]),'trial_data');
            if isfield(trial_data(1),'result')
                trial_data = trial_data(getTDidx(trial_data,'result',{'R'}));
            end
            
            if filter_trials
                [~,bad_idx] = filterTrials(trial_data);
                good_trial_idx = ~bad_idx;
            end
            
            metric = abs(get_learning_metrics(trial_data, 'angle')*180/pi);
            %metric = abs(get_learning_metrics(trial_data, 'corr'));%*180/pi);
            
            for e = 1:length(epochs)
                if file == 1
                    behavior{e} = repmat(metric(e_inds{e})',size(temp_metric_cv,1),1);
                    if filter_trials
                        good_trials{e} = repmat(good_trial_idx(e_inds{e}),size(temp_metric_cv,1),1);
                    end
                else
                    behavior{e} = [behavior{e}(:,1:length(e_inds{e})); repmat(metric(e_inds{e})',size(temp_metric_cv,1),1)];
                    if filter_trials
                        good_trials{e} = [good_trials{e}(:,1:length(e_inds{e})); repmat(good_trial_idx(e_inds{e}),size(temp_metric_cv,1),1)];
                    end
                end
            end
        end
    else
        disp([filenames{file} ' not found...']);
    end
end

out_struct.cv = cv;
out_struct.e_inds = e_inds;
out_struct.e_pr2 = e_pr2;
out_struct.good_cells = good_cells;
out_struct.total_significant = total_significant;
out_struct.total_cells = total_cells;

if do_behavior
    out_struct.behavior = behavior;
    if filter_trials
        out_struct.good_trials = good_trials;
    end
end

out_struct.params = params;