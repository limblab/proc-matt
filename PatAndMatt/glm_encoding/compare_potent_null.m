%%
clear
clc;
close all;

dataSummary;

basenames = {'trainad','trainad'};
extranames = {'potent_move','null_move'};
array_models = {'PMd-M1','PMd-M1'};

sessions = { ...
    'Chewie','2016-09-15'; ... % CF
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

pert = 'FF';
tasks = {'CO'};
dates = sessions(:,2);
monkeys = unique(sessions(:,1));

which_metric = 'rpr2'; % 'rpr2','pr2_full','pr2_basic'
pr2_cutoff = 0.01;
pr2_op = 'min'; % which operation for filtering ('min','max','mean','median')
pr2_wo_check = false; % only keeps cells that predict in WO

plot_op = 'mean';
group_size = 5;

do_norm = true;
do_diff = true;
error_bars = 'ste'; % 'boot','ste'
num_bootstraps = 2000;
how_to_group = 'average';

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
        filepaths{file} = fullfile(rootDir,TDDir,outputSubdir,[pert '-' array_models{idx_cond} '_' filenames{file} '.mat']);
    end
    
    out_struct = get_plot_metrics(filepaths, ...
        struct( ...
        'which_metric',which_metric, ...
        'epochs',{epochs}, ...
        'pr2_cutoff',pr2_cutoff, ...
        'pr2_op',pr2_op, ...
        'pr2_wo_check', pr2_wo_check, ...
        'do_good_cells',false, ...
        'do_behavior',true));
    
    cv = out_struct.cv;
    e_pr2 = out_struct.e_pr2;
    e_inds = out_struct.e_inds;
    good_cells = out_struct.good_cells;
    behavior = out_struct.behavior;
    
    for e = 1:length(epochs)
        v = e_pr2{e};
        b = behavior{e};
        if ~isempty(v)
            if do_diff
                v = v - repmat(mean(cv,2),1,size(v,2));
            end
            if do_norm
                v = v ./ repmat(abs(mean(cv,2)),1,size(v,2));
            end
            
            % group together some number of trials
            if group_size > 1
                group_idx = 1:group_size:size(v,2);
                
                switch lower(how_to_group)
                    case 'average'
                        temp_v = zeros(size(v,1),length(group_idx)-1);
                        temp_b = zeros(size(b,1),length(group_idx)-1);
                        for i = 1:length(group_idx)-1
                            temp_v(:,i) = mean(v(:,group_idx(i):group_idx(i+1)-1),2);
                            temp_b(:,i) = mean(b(:,group_idx(i):group_idx(i+1)-1),2);
                        end
                    case 'pool'
                        temp_v = zeros(group_size*size(v,1),length(group_idx)-1);
                        temp_b = zeros(group_size*size(b,1),length(group_idx)-1);
                        for i = 1:length(group_idx)-1
                            temp_v(:,i) = reshape(v(:,group_idx(i):group_idx(i+1)-1),group_size*size(v,1),1);
                            temp_b(:,i) = reshape(b(:,group_idx(i):group_idx(i+1)-1),group_size*size(b,1),1);
                        end
                end
                
                v = temp_v;
                b = temp_b;
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
good_idx = all_good_cells{1} & all_good_cells{2};

d = all_r2{2}(good_idx,:);% - all_r2{1}(good_idx,:);
err = all_err{1}(good_idx,:);

switch lower(plot_op)
    case 'median'
        m = median(d,1);
        switch lower(error_bars)
            case 'ste'
                s = [m-std(d,[],1)./sqrt(size(d,1)); m+std(d,[],1)./sqrt(size(d,1))];
            case 'boot'
                s = zeros(2,size(d,2));
                for k = 1:size(d,2)
                    temp = d(:,k);
                    bs = randi(length(temp),length(temp),num_bootstraps);
                    s(:,k) = prctile( median(temp(bs),1) ,[2.5,97.5])';
                end
        end
    case 'mean'
        m = mean(d,1);
        switch lower(error_bars)
            case 'ste'
                s = [m-std(d,[],1)./sqrt(size(d,1)); m+std(d,[],1)./sqrt(size(d,1))];
            case 'boot'
                s = zeros(2,size(d,2));
                for k = 1:size(d,2)
                    temp = d(:,k);
                    bs = randi(length(temp),length(temp),num_bootstraps);
                    s(:,k) = prctile( mean(temp(bs),1) ,[2.5,97.5])';
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
title([extranames{2} ' - ' extranames{1}],'FontSize',14);



%% Correlate difference with error
m_x = median(err,1);
[~,I] = sort(m_x);
s_x = [m_x - std(err,[],1)./sqrt(size(err,1)); m_x + std(err,[],1)./sqrt(size(err,1))];

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
xlabel('Average Error in Movements (Deg)','FontSize',14);
ylabel('Diff of Potent and Null p-R2','FontSize',14);
title([extranames{2} ' - ' extranames{1}],'FontSize',14);




