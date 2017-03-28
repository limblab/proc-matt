%%
clear
clc;
close all;

dataSummary;

basenames = {'trainad'};
extranames = {'singleneurons'};
array_models = {'PMd-M1'};

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
        };

pert = 'FF';
tasks = {'CO'};
dates = sessions(:,2);
monkeys = unique(sessions(:,1));

which_metric = 'rpr2'; % 'rpr2','pr2_full','pr2_basic'
pr2_cutoff = 0.0;
pr2_op = 'min'; % which operation for filtering ('min','max','mean','median')
pr2_ad_check = false; % only keeps cells that predict in WO
basic_pr2_check = false;

plot_op = 'mean';
group_size = 5;

do_norm = true;
do_diff = true;
error_bars = 'ste'; % 'boot','ste'
num_bootstraps = 1000;
how_to_group = 'pool';

epochs = {'AD'};


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
[all_r2, all_cv, all_err, all_good_cells] = deal(cell(1,length(basenames)));
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
        'do_good_cells',false, ...
        'do_behavior',true, ...
        'basic_pr2_check',basic_pr2_check));
    
    cv = out_struct.cv;
    e_pr2 = out_struct.e_pr2;
    e_inds = out_struct.e_inds;
    good_cells = out_struct.good_cells;
    behavior = out_struct.behavior;
    
    for e = 1:length(epochs)
        % get the total behavioral change
        %   abs(first 10 trials - last 10 trials)
        
        b = behavior{e};
        b_1 = abs(nanmean(b(:,1:5),2));
        b_2 = abs(nanmean(b(:,25:30),2));
        
        v = e_pr2{e}(:,1:41);
        b = behavior{e}(:,1:41);
        if ~isempty(v)
            if do_diff
                v = v - repmat(mean(cv,2),1,size(v,2));
            end
            if do_norm
                v = v ./ repmat(abs(mean(cv,2)),1,size(v,2));
            end
            
            bad_idx = abs(b) > 3 * std(reshape(b,numel(b),1)) | abs(v) > 3 * std(reshape(v,numel(v),1));
            if sum(bad_idx) > 0
                disp('outliers');
            end
            
            v(bad_idx) = NaN;
            b(bad_idx) = NaN;
            
            % group together some number of trials
            if group_size > 1
                group_idx = 1:group_size:size(v,2);
                
                switch lower(how_to_group)
                    case 'average'
                        temp_v = zeros(size(v,1),length(group_idx)-1);
                        temp_b = zeros(size(b,1),length(group_idx)-1);
                        for i = 1:length(group_idx)-1
                            temp_v(:,i) = nanmean(v(:,group_idx(i):group_idx(i+1)-1),2);
                            temp_b(:,i) = nanmean(b(:,group_idx(i):group_idx(i+1)-1),2);
                        end
                    case 'pool'
                        temp_v = zeros(group_size*size(v,1),length(group_idx)-1);
                        temp_b = zeros(group_size*size(b,1),length(group_idx)-1);
                        for i = 1:length(group_idx)-1
                            temp_v(:,i) = reshape(v(:,group_idx(i):group_idx(i+1)-1),group_size*size(v,1),1);
                            temp_b(:,i) = reshape(b(:,group_idx(i):group_idx(i+1)-1),group_size*size(b,1),1);
                        end
                        b_1 = repmat(b_1,group_size,1);
                        b_2 = repmat(b_2,group_size,1);
                    case 'slide'
                        group_idx = 1:size(v,2)-group_size;
                        temp_v = zeros(size(v,1),length(group_idx));
                        temp_b = zeros(size(b,1),length(group_idx));
                        for i = group_idx
                            temp_v(:,i) = nanmean(v(:,i:i+group_size),2);
                            temp_b(:,i) = nanmean(b(:,i:i+group_size),2);
                        end
                end
                
                v = temp_v;
                b = temp_b;%(repmat(b_1,1,size(temp_b,2)) - temp_b)./repmat(abs(b_1),1,size(temp_b,2));
            end
        end
    end
    all_r2{idx_cond} = v;
    all_cv{idx_cond} = cv;
    all_err{idx_cond} = b;
    all_good_cells{idx_cond} = good_cells;
end

%% Plot difference between null and potent
%filter out neurons
good_idx = all_good_cells{1}==1;% & all_good_cells{2};

d = all_r2{1}(good_idx,:);% - all_r2{1}(good_idx,:);
err = all_err{1}(good_idx,:);

switch lower(plot_op)
    case 'median'
        m = nanmedian(d,1);
        switch lower(error_bars)
            case 'ste'
                s = [m-nanstd(d,[],1)./sqrt(size(d,1)); m+nanstd(d,[],1)./sqrt(size(d,1))];
            case 'boot'
                s = zeros(2,size(d,2));
                for k = 1:size(d,2)
                    temp = d(:,k);
                    bs = randi(length(temp),length(temp),num_bootstraps);
                    s(:,k) = prctile( nanmedian(temp(bs),1) ,[2.5,97.5])';
                end
        end
    case 'mean'
        m = nanmean(d,1);
        switch lower(error_bars)
            case 'ste'
                s = [m-nanstd(d,[],1)./sqrt(size(d,1)); m+nanstd(d,[],1)./sqrt(size(d,1))];
            case 'boot'
                s = zeros(2,size(d,2));
                for k = 1:size(d,2)
                    temp = d(:,k);
                    bs = randi(length(temp),length(temp),num_bootstraps);
                    s(:,k) = prctile( nanmean(temp(bs),1) ,[2.5,97.5])';
                end
        end
end

all_x = group_size*(1:size(d,2)) - (group_size/2);

figure;
hold all;
plot(all_x,m,'LineWidth',2,'Color',plot_colors(1,:));
patch([all_x, fliplr(all_x)],[s(1,:), fliplr(s(2,:))],plot_colors(1,:),'FaceAlpha',0.3,'EdgeAlpha',0.5,'EdgeColor',plot_colors(1,:));
set(gca,'Box','off','TickDir','out','FontSize',14);
axis('tight');
xlabel('Trials from start of learning','FontSize',14);
ylabel('Difference in relative pseudo-r2','FontSize',14);



%% Correlate difference with error
m_x = zeros(1,size(err,2));
s_x = zeros(2,size(err,2));
for k = 1:size(err,2)
    m_x(k) = nanmean(err(:,k),1);
    s_x(:,k) = [m_x(k) - nanstd(err(:,k),[],1)./sqrt(size(err(:,k),1)); m_x(k) + nanstd(err(:,k),[],1)./sqrt(size(err(:,k),1))];
end

% m_x = nanmean(err,1);
%         switch lower(error_bars)
%             case 'ste'
%                 s_x = [m_x - nanstd(err,[],1)./sqrt(size(err,1)); m_x + nanstd(err,[],1)./sqrt(size(err,1))];
%             case 'boot'
%                 s_x = zeros(2,size(err,2));
%                 for k = 1:size(err,2)
%                     temp = unique(err(:,k));
%                     bs = randi(length(temp),length(temp),num_bootstraps);
%                     s_x(:,k) = prctile( nanmean(temp(bs),1) ,[2.5,97.5])';
%                 end
%         end
        
figure;
hold all;

for i = 1:length(m_x)
    plot(m_x(i),m(i),'o','LineWidth',2,'Color',plot_colors(1,:))
    plot([s_x(1,i), s_x(2,i)],[m(i), m(i)],'-','LineWidth',1,'Color',plot_colors(1,:));
    plot([m_x(i), m_x(i)],[s(1,i),s(2,i)],'-','LineWidth',1,'Color',plot_colors(1,:));
end
axis('square');
set(gca,'XLim',[floor(min(s_x(1,:))), ceil(max(s_x(2,:)))]);
V = axis;

% fit line to means
[b,~,~,~,stats] = regress(m',[ones(numel(m),1), m_x']);
plot(V(1):V(2),b(1)+b(2)*(V(1):V(2)),'k-','LineWidth',3);

text(0.75*V(2),0.7*V(4),['r = ' num2str(sqrt(stats(1)),2) '; p = ' num2str(stats(3),2)],'FontSize',12);

set(gca,'Box','off','TickDir','out','FontSize',14);
axis('tight');
axis('square');
% xlabel('% Adaptation','FontSize',14);
xlabel('Error (Deg)','FontSize',14);
ylabel('Change in rpR2','FontSize',14);




