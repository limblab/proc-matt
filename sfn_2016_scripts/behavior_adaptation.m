clear; clc; close all;
dataSummary;

sessions = { ...
%     'Chewie','2016-09-09'; ... % VR
%     'Chewie','2016-09-12'; ...
%     'Chewie','2016-09-14'; ...
%     'Chewie','2016-10-06'; ...
    'Mihili','2014-03-03'; ...
    'Mihili','2014-03-04'; ...
    'Mihili','2014-03-06'; ...
    'Mihili','2015-06-23'; ...
    'Mihili','2015-06-25'; ...
    'Mihili','2015-06-26'; ...
%     'Chewie','2016-09-15'; ... % CF
    'Chewie','2016-09-19'; ... % CF
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
%             'MrT','2013-08-19'; ... % CF
%             'MrT','2013-08-21'; ...
%             'MrT','2013-08-23'; ...
%             'MrT','2013-09-03'; ... %VR
%             'MrT','2013-09-05'; ...
%             'MrT','2013-09-09'; ...
    %     };
    % % CHEWIE M1 ONLY
    % % sessions = { ...
    'Chewie','2015-06-30'; ...
    'Chewie','2015-07-01'; ...
    'Chewie','2015-07-03'; ...
    'Chewie','2015-07-06'; ...
    'Chewie','2015-07-10'; ...
    'Chewie','2015-07-14'; ...
    'Chewie','2015-07-15'; ...
    'Chewie','2015-07-16'; ...
    };

pert = 'FF';
tasks = {'CO'};
dates = sessions(:,2);
monkeys = unique(sessions(:,1));

session_idx = find(ismember(filedb.Monkey,monkeys) & ismember(filedb.Perturbation,pert) & ismember(filedb.Task,tasks) & ismember(filedb.Date,dates));
filenames = cell(1,length(session_idx));
for s = 1:length(session_idx)
    filenames{s} = [filedb.Monkey{session_idx(s)} '_' filedb.Task{session_idx(s)} '_' filedb.Perturbation{session_idx(s)} '_' filedb.Date{session_idx(s)}];
end


epochs = {'BL','AD','WO'};

which_metric = 'angle';
remove_outliers = true;

%
switch lower(which_metric)
    case 'angle'
        err = cell(length(filenames),3);
        for file = 1:length(filenames)
            load(['F:\TrialDataFiles\' filenames{file} '.mat']);
            
            %     td = trial_data;
            td = filterTrials(trial_data);
            temp = get_learning_metrics(td,'angle');
            for i = 1:length(epochs)
                idx = getTDidx(td,'epoch',epochs{i});
                err{file,i} = sign(td(1).perturbation_info(end))*temp(idx);
            end
        end
    case 'corr'
        err = cell(length(filenames),3);
        for file = 1:length(filenames)
            load(['F:\TrialDataFiles\' filenames{file} '.mat']);
            
            %     td = trial_data;
            td = filterTrials(trial_data);
            temp = get_learning_metrics(td,'corr');
            for i = 1:length(epochs)
                idx = getTDidx(td,'epoch',epochs{i});
                err{file,i} = temp(idx);
            end
        end
end

%%
switch lower(which_metric)
    case 'angle'
        close all;
        % n_trials = 88;
        % n_trials = min(cellfun(@length,reshape(err,numel(err),1)))-group_size;
        group_size = 20;
        figure; subplot1(1,3);
        
        if remove_outliers
            for file = 1:size(err,1)
                temp = cat(1,err{file,:});
                for i = 1:size(err,2)
                    temp2 = err{file,i};
                    if  sum( abs(temp2) > 7*circular_std(reshape(temp,numel(temp),1))) > 0
                        disp('outliers');
                    end
                    temp2( abs(temp2) > 7*circular_std(reshape(temp,numel(temp),1))) = NaN;
                    err{file,i} = temp2;
                end
            end
        end
        
        for i = 1:length(epochs)
            subplot1(i);
            n_trials = min([min(cellfun(@length,err(:,i)))-group_size,500]);
            temp = zeros(length(filenames),n_trials);
            for j = 1:length(filenames)
                for k = 1:n_trials
                    temp(j,k) = circular_mean(err{j,i}(k:k+group_size));
                end
            end
            if strcmpi(epochs{i},'ad') && strcmpi(pert,'vr');
                temp = temp + abs(td(1).perturbation_info(end));
            end
            m = circular_mean(temp,[],1)*180/pi;
            if 0 % bootstrap
                s = zeros(2,size(temp,2));
                for k = 1:size(temp,2)
                    bs = randi(size(temp,1),size(temp,1),1000);
                    fuck = temp(:,k);
                    s(:,k) = 180/pi*prctile(circular_mean(fuck(bs),[],1),[2.5,97.5]);
                end
                patch([1:n_trials, fliplr(1:n_trials)],[s(1,:),fliplr(s(2,:))],[0.7 0.7 0.7]);
            else % std err
                s = 180/pi*std(temp,[],1)/sqrt(size(temp,1));
                patch([1:n_trials, fliplr(1:n_trials)],[m-s,fliplr(m+s)],[0.7 0.7 0.7]);
            end
            
            plot(m,'k','LineWidth',3);
            set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[-30 30],'XLim',[1,n_trials]);
        end
        
    case 'corr'
        close all;
        % n_trials = 88;
        % n_trials = min(cellfun(@length,reshape(err,numel(err),1)))-group_size;
        group_size = 10;
        figure; hold all;
        x_count = 0;
        
        if remove_outliers
            for file = 1:size(err,1)
                temp = cat(1,err{file,:});
                for i = 1:size(err,2)
                    temp2 = err{file,i};
                    if  sum(temp - mean(temp) > 2*std(reshape(temp,numel(temp),1))) > 0
                        disp('outliers');
                    end
                    temp2( temp2 - mean(temp) > 2*std(reshape(temp,numel(temp),1))) = NaN;
                    err{file,i} = temp2;
                end
            end
        end
        
        for i = 1:length(epochs)
            n_trials = min([min(cellfun(@length,err(:,i)))-group_size,70]);
            temp = zeros(length(filenames),n_trials);
            for j = 1:length(filenames)
                for k = 1:n_trials
                    temp(j,k) = mean(err{j,i}(k:k+group_size));
                end
            end
            
            m = nanmean(temp,1);
            if 1 % bootstrap
                s = zeros(2,size(temp,2));
                for k = 1:size(temp,2)
                    bs = randi(size(temp,1),size(temp,1),1000);
                    fuck = temp(:,k);
                    s(:,k) = prctile(nanmean(fuck(bs),1),[2.5,97.5]);
                end
                patch(x_count+[1:n_trials, fliplr(1:n_trials)],[s(1,:),fliplr(s(2,:))],[0.7 0.7 0.7]);
            else % std err
                s = nanstd(temp,[],1)/sqrt(size(temp,1));
                patch(x_count+[1:n_trials, fliplr(1:n_trials)],[m-s,fliplr(m+s)],[0.7 0.7 0.7]);
            end
            plot(x_count+(1:length(m)),m,'k','LineWidth',2);
            x_count = x_count + length(m);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[0.3 1]);
end


%%

clear; clc; close all;
dataSummary;

sessions = { ...
%     'Chewie','2016-09-09'; ... % VR
%     'Chewie','2016-09-12'; ...
%     'Chewie','2016-09-14'; ...
%     'Chewie','2016-10-06'; ...
    'Mihili','2014-03-03'; ...
    'Mihili','2014-03-04'; ...
    'Mihili','2014-03-06'; ...
    'Mihili','2015-06-23'; ...
    'Mihili','2015-06-25'; ...
    'Mihili','2015-06-26'; ...
            'MrT','2013-09-03'; ... %VR
            'MrT','2013-09-05'; ...
            'MrT','2013-09-09'; ...
    %     };
    % % CHEWIE M1 ONLY
    % % sessions = { ...
    'Chewie','2015-06-30'; ...
    'Chewie','2015-07-01'; ...
    'Chewie','2015-07-03'; ...
    'Chewie','2015-07-06'; ...
    'Chewie','2015-07-10'; ...
    'Chewie','2015-07-14'; ...
    'Chewie','2015-07-15'; ...
    'Chewie','2015-07-16'; ...
    };

pert = 'VR';
tasks = {'CO'};
dates = sessions(:,2);
monkeys = unique(sessions(:,1));

session_idx = find(ismember(filedb.Monkey,monkeys) & ismember(filedb.Perturbation,pert) & ismember(filedb.Task,tasks) & ismember(filedb.Date,dates));
filenames = cell(1,length(session_idx));
for s = 1:length(session_idx)
    filenames{s} = [filedb.Monkey{session_idx(s)} '_' filedb.Task{session_idx(s)} '_' filedb.Perturbation{session_idx(s)} '_' filedb.Date{session_idx(s)}];
end


epochs = {'BL','AD','WO'};

which_metric = 'angle';
remove_outliers = true;

switch lower(which_metric)
    case 'angle'
        err = cell(length(filenames),3);
        for file = 1:length(filenames)
            load(['F:\TrialDataFiles\' filenames{file} '.mat']);
            
            %     td = trial_data;
            td = filterTrials(trial_data);
            temp = get_learning_metrics(td,'angle');
            for i = 1:length(epochs)
                idx = getTDidx(td,'epoch',epochs{i});
                err{file,i} = sign(td(1).perturbation_info(end))*temp(idx);
            end
        end
    case 'corr'
        err = cell(length(filenames),3);
        for file = 1:length(filenames)
            load(['F:\TrialDataFiles\' filenames{file} '.mat']);
            
            %     td = trial_data;
            td = filterTrials(trial_data);
            temp = get_learning_metrics(td,'corr');
            for i = 1:length(epochs)
                idx = getTDidx(td,'epoch',epochs{i});
                err{file,i} = temp(idx);
            end
        end
end
%%
switch lower(which_metric)
    case 'angle'
        close all;
        % n_trials = 88;
        % n_trials = min(cellfun(@length,reshape(err,numel(err),1)))-group_size;
        group_size = 5;
        figure; subplot1(1,3);
        for i = 1:length(epochs)
            subplot1(i);
            n_trials = min([min(cellfun(@length,err(:,i)))-group_size,500]);
            temp = zeros(length(filenames),n_trials);
            for j = 1:length(filenames)
                for k = 1:n_trials
                    temp(j,k) = circular_mean(err{j,i}(k:k+group_size));
                end
            end
            if strcmpi(epochs{i},'ad') && strcmpi(pert,'vr');
                temp = temp + abs(td(1).perturbation_info(end));
            end
            m = circular_mean(temp,[],1)*180/pi;
            if 1 % bootstrap
                s = zeros(2,size(temp,2));
                for k = 1:size(temp,2)
                    bs = randi(size(temp,1),size(temp,1),1000);
                    fuck = temp(:,k);
                    s(:,k) = 180/pi*prctile(circular_mean(fuck(bs),[],1),[2.5,97.5]);
                end
                patch([1:n_trials, fliplr(1:n_trials)],[s(1,:),fliplr(s(2,:))],[0.7 0.7 0.7]);
            else % std err
                s = 180/pi*std(temp,[],1)/sqrt(size(temp,1));
                patch([1:n_trials, fliplr(1:n_trials)],[m-s,fliplr(m+s)],[0.7 0.7 0.7]);
            end
            
            plot(m,'k','LineWidth',3);
            set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[-30 30],'XLim',[1,n_trials]);
        end
        
    case 'corr'
        close all;
        % n_trials = 88;
        % n_trials = min(cellfun(@length,reshape(err,numel(err),1)))-group_size;
        group_size = 2;
        figure; hold all;
        x_count = 0;
        
        if remove_outliers
            for file = 1:size(err,1)
                temp = cat(1,err{file,:});
                for i = 1:size(err,2)
                    temp2 = err{file,i};
                    if  sum(temp - mean(temp) > 2*std(reshape(temp,numel(temp),1))) > 0
                        disp('outliers');
                    end
                    temp2( temp2 - mean(temp) > 2*std(reshape(temp,numel(temp),1))) = NaN;
                    err{file,i} = temp2;
                end
            end
        end
        
        for i = 1:length(epochs)
            n_trials = min([min(cellfun(@length,err(:,i)))-group_size,70]);
            temp = zeros(length(filenames),n_trials);
            for j = 1:length(filenames)
                for k = 1:n_trials
                    temp(j,k) = mean(err{j,i}(k:k+group_size));
                end
            end
            
            m = nanmean(temp,1);
            if 1 % bootstrap
                s = zeros(2,size(temp,2));
                for k = 1:size(temp,2)
                    bs = randi(size(temp,1),size(temp,1),1000);
                    fuck = temp(:,k);
                    s(:,k) = prctile(nanmean(fuck(bs),1),[2.5,97.5]);
                end
                patch(x_count+[1:n_trials, fliplr(1:n_trials)],[s(1,:),fliplr(s(2,:))],[0.7 0.7 0.7]);
            else % std err
                s = nanstd(temp,[],1)/sqrt(size(temp,1));
                patch(x_count+[1:n_trials, fliplr(1:n_trials)],[m-s,fliplr(m+s)],[0.7 0.7 0.7]);
            end
            plot(x_count+(1:length(m)),m,'k','LineWidth',2);
            x_count = x_count + length(m);
        end
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[0.3 1]);
end












