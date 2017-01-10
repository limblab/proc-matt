%%
clear
% clc;
% close all;

if 0
    a = reshape(e_pr2{2},numel(e_pr2{2}),1);
    g = reshape( repmat(1:size(e_pr2{2},2),size(e_pr2{2},1),1), numel(e_pr2{2}), 1 );
    anovan(a,g,'Display','off')
end

dataSummary;

outputSubdir = 'trainad_null8';

% sessions = { ...
%     'Chewie','2016-09-09'; ... % VR
%     'Chewie','2016-09-12'; ...
%     'Chewie','2016-09-14'; ...
%     'Mihili','2014-03-03'; ...
%     'Mihili','2014-03-04'; ...
%     'Mihili','2014-03-06'; ...
%     'Mihili','2015-06-23'; ...
%     'Mihili','2015-06-25'; ...
%     'Mihili','2015-06-26'; ...
%     'Chewie','2016-09-15'; ... % CF
%     'Chewie','2016-09-19'; ...
%     'Chewie','2016-10-05'; ...
%     'Chewie','2016-10-07'; ...
%     'Mihili','2014-02-03'; ...
%     'Mihili','2014-02-17'; ...
%     'Mihili','2014-02-18'; ...
%     'Mihili','2014-03-07'; ...
%     'Mihili','2015-06-15'; ...
%     'Mihili','2015-06-16'; ...
%     'Mihili','2015-06-17'; ...
%     'MrT','2013-08-19'; ... % CF
%     'MrT','2013-08-21'; ...
%     'MrT','2013-08-23'; ...
%     'MrT','2013-09-03'; ... % VR
%     'MrT','2013-09-05'; ...
%     'MrT','2013-09-09'; ...
%     };

sessions = { ...
    'Chewie','2016-09-15'; ... % CF
    'Chewie','2016-09-19'; ...
    'Chewie','2016-10-05'; ...
    'Chewie','2016-10-07'; ...
    'Mihili','2014-02-03'; ...
    'Mihili','2014-02-17'; ...
    'Mihili','2014-02-18'; ...
    'Mihili','2014-03-07'; ...
%     'Chewie','2016-09-09'; ... % VR
%     'Chewie','2016-09-12'; ...
%     'Chewie','2016-09-14'; ...
%     'Mihili','2014-03-03'; ...
%     'Mihili','2014-03-04'; ...
%     'Mihili','2014-03-06'; ...
};


perts = {'FF'};
tasks = {'CO'};
dates = sessions(:,2);
monkeys = unique(sessions(:,1));

% array_pairs = {'M1-M1','PMd-PMd','PMd-M1'};
array_pairs = {'PMd-M1'};

which_metric = 'rpr2'; % 'rpr2','pr2_full','pr2_basic'
pr2_cutoff = 0.01;
pr2_op = 'mean'; % which operation for filtering ('min','max','mean','median')
pr2_wo_check = false; % only keeps cells that predict in WO

plot_op = 'median';
group_size = 10;

do_norm = true;
do_diff = true;
error_bars = 'boot'; % 'boot','ste'
num_bootstraps = 1500;

do_regression_line = false;
add_plot = true;
do_subplot = false;

epochs = {'BL','AD','WO'};

if do_norm
    min_y = -2; max_y = 0.5;
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

session_idx = ismember(filedb.Monkey,monkeys) & ismember(filedb.Perturbation,perts) & ismember(filedb.Task,tasks) & ismember(filedb.Date,dates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% max_y = 0;
% min_y = 0;
min_vaf = 0;
max_vaf = 0;
all_shit = [];
for idx_pert = 1:length(perts)
    if ~add_plot
        figure;
    end
    idx = find(strcmpi(perts{idx_pert},filedb.Perturbation) & session_idx);
    filenames = cell(1,length(idx));
    for s = 1:length(idx)
        filenames{s} = [filedb.Monkey{idx(s)} '_' filedb.Task{idx(s)} '_' filedb.Perturbation{idx(s)} '_' filedb.Date{idx(s)}];
    end
    
    for model = 1:length(array_pairs)
        if do_subplot
            subplot(1,length(array_pairs),model); hold all;
        else
            hold all;
        end
        
        % build list of filenames
        filepaths = cell(1,length(filenames));
        for file = 1:length(filenames)
            filepaths{file} = fullfile(rootDir,TDDir,outputSubdir,[perts{idx_pert} '-' array_pairs{model} '_' filenames{file} '.mat']);
        end
        
        out_struct = get_plot_metrics(filepaths,struct('which_metric',which_metric,'epochs',{epochs},'pr2_cutoff',pr2_cutoff,'pr2_op',pr2_op,'pr2_wo_check',pr2_wo_check));
        cv = out_struct.cv;
        e_pr2 = out_struct.e_pr2;
        e_inds = out_struct.e_inds;
        total_significant = out_struct.total_significant;
        total_cells = out_struct.total_cells;
        
        disp([array_pairs{model} ' - % cells with significant rel-pseudo-r2: ' num2str(total_significant) '/' num2str(total_cells)]);
        disp([array_pairs{model} ' - Mean baseline metric: ' num2str(mean(mean(cv,2))) ' +/- ' num2str(std(mean(cv,2)))]);
        
        [all_x, all_y, all_y_std] = deal([]);
        for e = 1:length(epochs)
            v = e_pr2{e};
            if ~isempty(v)
                if do_diff
                    v = v - repmat(mean(cv,2),1,size(v,2));
                    if do_norm
                        v = v ./ repmat(abs(mean(cv,2)),1,size(v,2));
                    end
                end
                
                % group together some number of trials
                if group_size > 1
                    group_idx = 1:group_size:size(v,2);
                    temp_v = zeros(group_size*size(v,1),length(group_idx)-1);
                    for i = 1:length(group_idx)-1
                        temp_v(:,i) = reshape(v(:,group_idx(i):group_idx(i+1)-1),group_size*size(v,1),1);
                    end
                    v = temp_v;
                end
                
                m = zeros(1,size(v,2));
                s = zeros(2,size(v,2));
                for k = 1:size(v,2)
                    switch lower(plot_op)
                        case 'median'
                            m(k) = median(v(:,k),1);
                            switch lower(error_bars)
                                case 'boot' % bootstrapped 95% conf int
                                    temp = v(:,k);
                                    bs = randi(length(temp),length(temp),num_bootstraps);
                                    s(:,k) = prctile( median(temp(bs),1) ,[2.5,97.5])';
                                case 'ste' % standard error of mean
                                    s(:,k) = [m(k)-std(v(:,k),1)./sqrt(size(v(:,k),1)); ...
                                        m(k)+std(v(:,k),1)./sqrt(size(v(:,k),1))];
                            end
                        case 'mean'
                            m(k) = mean(v(:,k),1);
                            switch lower(error_bars)
                                case 'boot' % bootstrapped 95% conf int
                                    temp = v(:,k);
                                    bs = randi(length(temp),length(temp),num_bootstraps);
                                    s(:,k) = prctile( mean(temp(bs),1) ,[2.5,97.5])';
                                case 'ste' % standard error of mean
                                    s(:,k) = [m(k)-std(v(:,k),1)./sqrt(size(v(:,k),1)); ...
                                        m(k)+std(v(:,k),1)./sqrt(size(v(:,k),1))];
                            end
                    end
                end
                
                if do_regression_line
                    %[b,~,~,~,stats] = regress(reshape(v,numel(v),1),[ones(numel(v),1), reshape(repmat(1:size(v,2),size(v,1),1),numel(v),1)]);
                    [b,~,~,~,stats] = regress(m',[ones(numel(m),1), (1:length(m))']);
                    plot(1+size(all_x,2)+(1:length(m)),b(1) + b(2)*(1:length(m)),'k-','LineWidth',2);
                    text(size(all_x,2),0,[num2str(stats(1),'%10.2e') ', ' num2str(stats(3),'%10.2e')]);
                end
                
                all_x = [all_x, 1+size(all_x,2)+(1:length(m))];
                all_y = [all_y, m];
                all_y_std = [all_y_std, s];
            end
        end
        
        all_x = all_x - 1;
        
        % min_y = min([min_y, all_y_std(1,:)]);
        % max_y = max([max_y, all_y_std(2,:)]);
        
        hold all;
        patch([all_x, fliplr(all_x)],[all_y_std(1,:), fliplr(all_y_std(2,:))],plot_colors(model,:),'FaceAlpha',0.3,'EdgeAlpha',0.5,'EdgeColor',plot_colors(model,:));
        h(model) = plot(all_x,all_y,'+','LineWidth',3,'Color',plot_colors(model,:));
        
    end
    
    % some end-of-script tweaking of stuff
    for model = 1:length(array_pairs)
        if do_subplot
            h = subplot(1,3,model);
            title(h,array_pairs{model});
        else
            h = gca;
        end
        axis('tight');
        set(h,'Box','off','TickDir','out','FontSize',14,'YLim',[min_y, max_y],'FontSize',14);
        %             'XTick',[floor(0+length(ad_inds)/2), ceil(1+length(ad_inds)+length(wo_inds)/2)],'XTickLabel',{'Force','Wash'},'FontSize',14);
        for e = 1:length(epochs)-1
            plot([sum(cellfun(@(x) floor(length(x)/group_size),e_inds(1:e))), sum(cellfun(@(x) floor(length(x)/group_size),e_inds(1:e)))]+0.5,[min_y,max_y],'k--','LineWidth',2);
        end
    end
    title(outputSubdir,'FontSize',14);
    if idx_pert == 1
        if do_norm
            ylabel(['norm change in ' which_metric],'FontSize',14);
        else
            ylabel(['change in ' which_metric],'FontSize',14);
        end
    end
    if idx_pert == length(perts)
        legend(h,array_pairs,'Box','off','FontSize',14,'Location','NorthEast');
    end
    
end