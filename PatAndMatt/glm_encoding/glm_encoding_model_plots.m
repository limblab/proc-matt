%%
clear
clc;
close all;

dataSummary;

outputSubdir = 'all';
sessions = { ...
    'Chewie','2016-09-15'; ...
%     'Chewie','2016-09-19'; ...
%     'Chewie','2016-10-05'; ...
%     'Chewie','2016-10-07'; ...
%     'Mihili','2014-02-17'; ...
%     'Mihili','2014-02-18'; ...
%     'Mihili','2014-03-07'; ...
%     'Mihili','2015-06-15'; ...
%     'Mihili','2015-06-16'; ...
%     'Mihili','2015-06-17'; ...
    };

% array_pairs = {'M1-M1','PMd-PMd','PMd-M1'};
array_pairs = {'PMd-PMd'};

which_metric = 'rpr2';
pr2_cutoff = 0;
pr2_op = 'min';
pr2_wo_check = true; % only keeps cells that predict in WO

do_norm = false;
do_diff = true;
do_conf_int = true; % bootstrapped 95% confidence bounds for error bars
num_bootstraps = 1000;

epochs = {'BL','AD','WO'};

if do_norm
    min_y = -5; max_y = 0.25;
    xmin = -5;
    xmax = 5;
    dx = 0.2;
else
    min_y = -0.75; max_y = 0.1;
    xmin = -0.3;
    xmax = 0.15;
    dx = 0.01;
end


%%
nbins = xmin-dx:dx:xmax+dx;

plot_colors = [0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    %     0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840];

session_idx = ismember(filedb.Monkey,sessions(:,1)) & ismember(filedb.Date,sessions(:,2));
perts = unique(filedb.Perturbation(session_idx));

figure;
if length(perts) > 1
    subplot1(1,length(perts));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% max_y = 0;
% min_y = 0;
min_vaf = 0;
max_vaf = 0;
all_shit = [];
for idx_pert = 1:length(perts)
    if length(perts) > 1
        subplot1(idx_pert);
    end
    
    idx = find(strcmpi(perts{idx_pert},filedb.Perturbation) & session_idx);
    filenames = cell(1,length(idx));
    for s = 1:length(idx)
        filenames{s} = [filedb.Monkey{idx(s)} '_' filedb.Task{idx(s)} '_' filedb.Perturbation{idx(s)} '_' filedb.Date{idx(s)}];
    end
    
    for model = 1:length(array_pairs)
        
        cv = [];
        % loop along files and group together
        [total_cells,total_significant] = deal(0);
        
        for file = 1:length(filenames) % loop along sessions
            load(fullfile(rootDir,TDDir,outputSubdir,[perts{idx_pert} '-' array_pairs{model} '_' filenames{file} '.mat']),'results','params');
            
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
            good_idx = temp_pr2(:,1) > pr2_cutoff;
            
            temp_metric = results.(which_metric);
            
            if pr2_wo_check % only take cells that can predict in the washout
                good_idx = good_idx & mean(mean(temp_metric(:,wo_inds(1:end),:),3),2) > pr2_cutoff;
            end
            
            temp_metric = temp_metric(good_idx,:,:);
            
            % get some stats
            total_significant = total_significant+sum(good_idx);
            total_cells = total_cells+length(good_idx);
            
            temp_metric_cv = results.([which_metric '_cv']);
            
            % parse out relevant params
            bl_inds = find(strcmpi(params.test_epochs,'BL'));
            ad_inds = find(strcmpi(params.test_epochs,'AD'));
            wo_inds = find(strcmpi(params.test_epochs,'WO'));
            
            if file > 1
                bl_inds = bl_inds(1:min([length(bl_inds), size(bl,2)]));
                ad_inds = ad_inds(1:min([length(ad_inds), size(ad,2)]));
                wo_inds = wo_inds(1:min([length(wo_inds), size(wo,2)]));
                bl = [bl(:,1:length(bl_inds)); mean(temp_metric(:,bl_inds,:),3)];
                ad = [ad(:,1:length(ad_inds)); mean(temp_metric(:,ad_inds,:),3)];
                wo = [wo(:,1:length(wo_inds)); mean(temp_metric(:,wo_inds,:),3)];
            else
                bl = mean(temp_metric(:,bl_inds,:),3);
                ad = mean(temp_metric(:,ad_inds,:),3);
                wo = mean(temp_metric(:,wo_inds,:),3);
            end
            
            cv = [cv; mean(temp_metric_cv(good_idx,:,:),3)];
        end
        
        disp([array_pairs{model} ' - % cells with significant rel-pseudo-r2: ' num2str(total_significant) '/' num2str(total_cells)]);
        
        disp([array_pairs{model} ' - Mean baseline metric: ' num2str(mean(mean(cv,2))) ' +/- ' num2str(std(mean(cv,2)))]);
        
        % OVER TIME
        plot_inds = [bl_inds;ad_inds;wo_inds];
        
        if do_diff
            bl = bl - repmat(mean(cv,2),1,size(bl,2));
            ad = ad - repmat(mean(cv,2),1,size(ad,2));
            wo = wo - repmat(mean(cv,2),1,size(wo,2));
        end
        
        if do_norm
            bl = bl ./ repmat(abs(mean(cv,2)),1,size(bl,2));
            ad = ad ./ repmat(abs(mean(cv,2)),1,size(ad,2));
            wo = wo ./ repmat(abs(mean(cv,2)),1,size(wo,2));
        end
        
        % BASELINE
        m = zeros(1,size(bl,2));
        s = zeros(2,size(bl,2));
        for k = 1:size(bl,2)
            m(k) = nanmedian(bl(:,k),1);
            if do_conf_int
                temp = bl(:,k);
                bs = zeros(1,num_bootstraps);
                for z = 1:num_bootstraps
                    bs(z) = nanmedian(temp(randi(length(temp),length(temp),1)),1);
                end
                s(:,k) = prctile(bs,[2.5,97.5])';
            else
                s(:,k) = [m(k)-nanstd(bl(:,k),1)./sqrt(size(bl(:,k),1)); ...
                    m(k)+nanstd(bl(:,k),1)./sqrt(size(bl(:,k),1))];
            end
        end
        all_x = 1:length(m);
        all_y = m;
        all_y_std = s;
        
        % ADAPTATION
        m = zeros(1,size(ad,2));
        s = zeros(2,size(ad,2));
        for k = 1:size(ad,2)
            m(k) = nanmedian(ad(:,k),1);
            if do_conf_int
                temp = ad(:,k);
                bs = zeros(1,num_bootstraps);
                for z = 1:num_bootstraps
                    bs(z) = nanmedian(temp(randi(length(temp),length(temp),1)),1);
                end
                s(:,k) = prctile(bs,[2.5,97.5])';
            else
                s(:,k) = [m(k)-nanstd(ad(:,k),1)./sqrt(size(ad(:,k),1)); ...
                    m(k)+nanstd(ad(:,k),1)./sqrt(size(ad(:,k),1))];
            end
        end
        all_x = [all_x, 1+size(bl,2)+(1:length(m))];
        all_y = [all_y, m];
        all_y_std = [all_y_std, s];
        
        % WASHOUT
        m = zeros(1,size(wo,2));
        s = zeros(2,size(wo,2));
        for k = 1:size(wo,2)
            m(k) = nanmedian(wo(:,k),1);
            if do_conf_int
                temp = wo(:,k);
                bs = zeros(1,num_bootstraps);
                for z = 1:num_bootstraps
                    bs(z) = nanmedian(temp(randi(length(temp),length(temp),1)),1);
                end
                s(:,k) = prctile(bs,[2.5,97.5])';
            else
                s(:,k) = [m(k)-nanstd(wo(:,k),1)./sqrt(size(wo(:,k),1)); ...
                    m(k)+nanstd(wo(:,k),1)./sqrt(size(wo(:,k),1))];
            end
        end
        all_x = [all_x, 1+size(bl,2)+size(ad,2)+(1:length(m))];
        all_y = [all_y, m];
        all_y_std = [all_y_std, s];
        
        all_x = all_x - 1;
        
        %         min_y = min([min_y, all_y_std(1,:)]);
        %         max_y = max([max_y, all_y_std(2,:)]);
        
        hold all;
        patch([all_x, fliplr(all_x)],[all_y_std(1,:), fliplr(all_y_std(2,:))],plot_colors(model,:),'FaceAlpha',0.3,'EdgeAlpha',0.5,'EdgeColor',plot_colors(model,:));
        h(model) = plot(all_x,all_y,'+','LineWidth',3,'Color',plot_colors(model,:));
        %         h(model) = plot(all_x,all_y,'-','LineWidth',3,'Color',plot_colors(model,:));
        
    end
    
    % some end-of-script tweaking of stuff
    for model = 1:length(array_pairs)
        set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 length(bl_inds)+length(ad_inds)+length(wo_inds)],'YLim',[min_y, max_y], ...
            'XTick',[floor(0+length(ad_inds)/2), ceil(1+length(ad_inds)+length(wo_inds)/2)],'XTickLabel',{'Force','Wash'},'FontSize',14);
        plot([length(bl_inds), length(bl_inds)]+0.5,[min_y, max_y],'k--','LineWidth',2);
        plot(0.5+length(bl_inds)+[length(ad_inds), length(ad_inds)],[min_y, max_y],'k--','LineWidth',2);
    end
    if idx_pert == 1
        if do_norm
            ylabel('Norm Change in Pseudo-R^2','FontSize',14);
        else
            ylabel('Change in Pseudo-R^2','FontSize',14);
        end
    end
    if idx_pert == length(perts)
        legend(h,array_pairs,'Box','off','FontSize',14,'Location','NorthEast');
    end
    
end