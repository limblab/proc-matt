%%
clear
clc;
close all;

dataSummary;

fuck_this = 9;
switch fuck_this
    case 1
        output_subdir = '10msec_bins';
        bl_inds = [1:2];
        ad_inds = [3:22];
        wo_inds = [23:42];
        which_files = 'all';
    case 2
        output_subdir = '20msec_bins';
        bl_inds = [1:2];
        ad_inds = [3:22];
        wo_inds = [23:42];
        which_files = 'all';
    case 3
        output_subdir = '50msec_bins_fine';
        bl_inds = [1];
        ad_inds = [2:21];
        wo_inds = [22:31];
        which_files = 'all';
    case 4
        output_subdir = 'new_standard_10bins';
        bl_inds = 1;
        ad_inds = [2:11];
        wo_inds = [12:21];
        which_files = 'all';
    case 5
        output_subdir = 'cv_epochs_force';
        bl_inds = 1;
        ad_inds = 2:3;
        wo_inds = 4:5;
        which_files = 'all';
    case 6
        output_subdir = 'allhistory_50msec_bins';
        bl_inds = [1];
        ad_inds = [2:21];
        wo_inds = [22:41];
        which_files = 'all';
    case 7
        output_subdir = '30msec_bins_fine';
        bl_inds = [1];
        ad_inds = [2:51];
        wo_inds = [52:101];
        which_files = 'all';
    case 8
        output_subdir = '30msec_bins';
        bl_inds = 1;
        ad_inds = [2:11];
        wo_inds = [12:21];
        which_files = 'all';
    case 9
        output_subdir = '50msec_bins_test';
        bl_inds = [];
        ad_inds = [1:10];
        wo_inds = [11:18];
        which_files = 'all';
    case 10
        output_subdir = 'single_trial';
        bl_inds = [];
        which_files = 2;
        switch which_files
            case 1
                ad_inds = [1:264];
                wo_inds = [265:523];
            case 2
                ad_inds = [1:201];
                wo_inds = [202:355];
        end
end

% 1: scatter plot
% 2: histogram difference for ff or vr, comparing washout and ad
% 3: over time
% 4: same as 2 but compares washout and ad
do_plot = 3;
pert = {'ff'}; % for 1 and 2

% array_pairs = {'M1-M1','PMd-PMd','PMd-M1'};
array_pairs = {'PMd-M1'};

min_rpr2 = 0;
min_pr2_full = 0;
which_metric = 'rpr2';

do_norm = false;
remove_outliers =true; % supported for 5 and 7 right now
outlier_alpha = 5;
do_conf_int = true; % bootstrapped 95% confidence bounds for error bars
num_bootstraps = 50;

use_cv = true; % if true, use cross validated baseline for diff
% note: if false, uses average of bl_inds entries in epoch predictions

if do_norm
    xmin = -5;
    xmax = 5;
    dx = 0.2;
else
    xmin = -0.15;
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

if ~iscell(pert)
    pert = {pert};
end

figure;
if do_plot ~= 3
    if length(pert) > 1
        error('Pick a perturbation dammit!');
    else
        subplot1(1,length(array_pairs));
    end
else
    if length(pert) > 1
        subplot1(1,length(pert));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_y = 0;
min_y = 0;
min_vaf = 0;
max_vaf = 0;
all_shit = [];
for idx_pert = 1:length(pert)
    if length(pert) > 1
        subplot1(idx_pert);
    end
    for i = 1:length(array_pairs)
        load(fullfile(rootDir,TDDir,output_subdir,[pert{idx_pert} '-' array_pairs{i}]));
        
%         which_files = strcmpi(params.use_files(:,4),'co');
        
        if ~ischar(which_files)
            results = results(which_files);
        end
        
        [all_m, all_m_bl,all_m_cv] = deal([]);
        % loop along files and group together
        [total_cells,total_significant] = deal(0);
        total_cells_session = zeros(1,length(results));
        for j = 1:length(results) % loop along sessions
            
            % get behavior adaptation
            
            if use_cv
                % temp_pr2 = max(results(j).pr2_full_cv(:,:,1),[],2);
                % temp_rpr2 = max(results(j).rpr2_cv(:,:,1),[],2);
                temp_pr2 = mean(mean(results(j).pr2_full_cv,3),2);
                temp_rpr2 = mean(mean(results(j).rpr2_cv,3),2);
                %good_idx = max(temp_pr2(:,:,1),[],2) > min_pr2_full & max(mean(temp_rpr2,3),[],2) > min_rpr2;
            else
                temp_pr2 = squeeze(results(j).pr2_full(:,bl_inds,:));
                temp_rpr2 = squeeze(results(j).rpr2(:,bl_inds,:));
            end
            
            good_idx = temp_pr2(:,1) > min_pr2_full & temp_rpr2(:,1) > min_rpr2;
            
            temp_metric = results(j).(which_metric);
            temp_metric = temp_metric(good_idx,:,:);
            
            % get some stats
            total_significant = total_significant+sum(good_idx);
            total_cells = total_cells+length(good_idx);
            total_cells_session(j) = length(good_idx);
            
            if ~isempty(bl_inds)
                temp_metric_bl = squeeze(mean(results(j).(which_metric)(:,bl_inds,:),2));
                temp_metric_bl = temp_metric_bl(good_idx,:,:);
            end
            
            
            if use_cv
                temp_metric_cv = results(j).([which_metric '_cv']);
                all_m_cv = [all_m_cv; mean(temp_metric_cv(good_idx,:,:),3)];
            else
                all_m_bl = [all_m_bl; mean(temp_metric_bl,3)];
            end
            % pool together metrics from all files for plotting
            all_m = [all_m; mean(temp_metric,3)];
        end
        
        disp([array_pairs{i} ' - % cells with significant rel-pseudo-r2: ' num2str(total_significant) '/' num2str(total_cells)]);
        
        if use_cv
            disp([array_pairs{i} ' - Mean baseline metric: ' num2str(mean(mean(all_m_cv,2))) ' +/- ' num2str(std(mean(all_m_cv,2)))]);
        else
            disp([array_pairs{i} ' - Mean baseline metric: ' num2str(mean(mean(all_m_bl,2))) ' +/- ' num2str(std(mean(all_m_bl,2)))]);
        end
        
        all_shit = [all_shit; mean(all_m_bl,2)];
        
        switch do_plot
            %%
            case 1
                min_vaf = min([min(all_m(:,[ad_inds,wo_inds]),[],1),0]);
                max_vaf = max(max(all_m(:,[ad_inds,wo_inds]),[],1));
                
                subplot(2,length(array_pairs),i);
                hold all;
                plot([min_vaf, max_vaf],[min_vaf,max_vaf],'k--','LineWidth',1);
                plot(mean(all_m_bl,2),all_m(:,ad_inds(1)),'o','LineWidth',1);
                plot(mean(all_m_bl,2),all_m(:,ad_inds(2)),'o','LineWidth',1);
                axis('square');
                set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[min_vaf max_vaf],'YLim',[min_vaf max_vaf]);
                if i==1
                    ylabel(['Adaptation ' which_metric],'FontSize',14);
                end
                title(array_pairs{i},'FontSize',16);
                
                subplot(2,length(array_pairs),i+length(array_pairs));
                hold all;
                plot([min_vaf, max_vaf],[min_vaf,max_vaf],'k--','LineWidth',1);
                plot(mean(all_m_bl,2),all_m(:,wo_inds(1)),'o','LineWidth',1);
                plot(mean(all_m_bl,2),all_m(:,wo_inds(2)),'o','LineWidth',1);
                axis('square');
                set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[min_vaf max_vaf],'YLim',[min_vaf max_vaf]);
                xlabel(['Baseline ' which_metric],'FontSize',14);
                if i == 1
                    ylabel(['Washout ' which_metric],'FontSize',14);
                end
                
            case 2
                %% HISTOGRAM DIFFERENCE VERSION
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                d_vaf = all_m - repmat(mean(all_m_bl,2),1,size(all_m,2));
                
                if do_norm
                    if use_cv
                        d_vaf = d_vaf ./ repmat(abs(mean(all_m_cv,2)),1,size(all_m,2));
                    else
                        d_vaf = d_vaf ./ repmat(abs(mean(all_m_bl,2)),1,size(all_m,2));
                    end
                end
                
                subplot(2,length(array_pairs),i);
                hold all;
                [N,X] = hist(d_vaf(:,ad_inds(1)),nbins);
                N = N./sum(N);
                bar(X,N,1,'FaceColor','r','FaceAlpha',0.7);
                [N,X] = hist(d_vaf(:,ad_inds(2)),nbins);
                N = N./sum(N);
                bar(X,N,1,'FaceAlpha',0.7);
                axis('tight'); V=axis;
                plot([0 0],V(3:4),'k--','LineWidth',2);
                set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[xmin xmax]);
                if i==1
                    switch lower(pert{idx_pert})
                        case 'ff'
                            ylabel(['Force ' which_metric],'FontSize',14);
                        case 'vr'
                            ylabel(['Rotation ' which_metric],'FontSize',14);
                    end
                end
                title(array_pairs{i},'FontSize',16);
                max_y = max([max_y, V(4)]);
                
                subplot(2,length(array_pairs),i+length(array_pairs));
                hold all;
                [N,X] = hist(d_vaf(:,wo_inds(1)),nbins);
                N = N./sum(N);
                bar(X,N,1,'FaceColor','r','FaceAlpha',0.7);
                [N,X] = hist(d_vaf(:,wo_inds(2)),nbins);
                N = N./sum(N);
                bar(X,N,1,'FaceAlpha',0.7);
                axis('tight'); V=axis;
                
                set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[xmin xmax]);
                xlabel('Baseline VAF','FontSize',14);
                if i == 1
                    ylabel(['Washout ' which_metric],'FontSize',14);
                end
                max_y = max([max_y, V(4)]);
                
                % do some statistical tests
                [~,p] = ttest2(d_vaf(:,ad_inds(2)),d_vaf(:,wo_inds(2)));
                disp([pert{idx_pert} ' ' array_pairs{i} ' - AD on WO: ' num2str(p)]);
                
            case 3
                %% OVER TIME
                plot_inds = [bl_inds,ad_inds,wo_inds];
                
                if do_norm
                    if use_cv
                        d_vaf = all_m - repmat(mean(all_m_cv,2),1,size(all_m,2));
                        d_vaf = d_vaf ./ repmat(abs(mean(all_m_cv,2)),1,size(all_m,2));
                    else
                        d_vaf = all_m - repmat(mean(all_m_bl,2),1,size(all_m,2));
                        d_vaf = d_vaf ./ repmat(abs(mean(all_m_bl,2)),1,size(all_m,2));
                    end
                else
                    if use_cv
                        d_vaf = all_m - repmat(mean(all_m_cv,2),1,size(all_m,2));
                    else
                        d_vaf = all_m - repmat(mean(all_m_bl,2),1,size(all_m,2));
                    end
                end
                
                bad_inds = zeros(size(d_vaf));
                if remove_outliers
                    for k = 1:size(d_vaf,2)
                        temp = prctile(d_vaf(:,k),[outlier_alpha/2, 100-outlier_alpha/2]);
                        bad_inds(:,k) = d_vaf(:,k) < temp(1) | d_vaf(:,k) > temp(2);
                    end
                end
                
                
                bl = d_vaf(:,bl_inds);
                ad = d_vaf(:,ad_inds);
                wo = d_vaf(:,wo_inds);
                
                % BASELINE
                m = zeros(1,size(bl,2));
                s = zeros(2,size(bl,2));
                for k = 1:size(bl,2)
                    m(k) = nanmedian(bl(~bad_inds(:,bl_inds(k)),k),1);
                    if do_conf_int
                        temp = bl(~bad_inds(:,bl_inds(k)),k);
                        bs = zeros(1,num_bootstraps);
                        for z = 1:num_bootstraps
                            bs(z) = nanmedian(temp(randi(length(temp),length(temp),1)),1);
                        end
                        s(:,k) = prctile(bs,[2.5,97.5])';
                    else
                        s(:,k) = [m(k)-nanstd(bl(~bad_inds(:,bl_inds(k)),k),1)./sqrt(size(bl(~bad_inds(:,bl_inds(k)),k),1)); ...
                            m(k)+nanstd(bl(~bad_inds(:,bl_inds(k)),k),1)./sqrt(size(bl(~bad_inds(:,bl_inds(k)),k),1))];
                    end
                end
                all_x = 1:length(m);
                all_y = m;
                all_y_std = s;
                
                % ADAPTATION
                m = zeros(1,size(ad,2));
                s = zeros(2,size(ad,2));
                for k = 1:size(ad,2)
                    m(k) = nanmedian(ad(~bad_inds(:,ad_inds(k)),k),1);
                    if do_conf_int
                        temp = ad(~bad_inds(:,ad_inds(k)),k);
                        bs = zeros(1,num_bootstraps);
                        for z = 1:num_bootstraps
                            bs(z) = nanmedian(temp(randi(length(temp),length(temp),1)),1);
                        end
                        s(:,k) = prctile(bs,[2.5,97.5])';
                    else
                        s(:,k) = [m(k)-nanstd(ad(~bad_inds(:,ad_inds(k)),k),1)./sqrt(size(ad(~bad_inds(:,ad_inds(k)),k),1)); ...
                            m(k)+nanstd(ad(~bad_inds(:,ad_inds(k)),k),1)./sqrt(size(ad(~bad_inds(:,ad_inds(k)),k),1))];
                    end
                end
                all_x = [all_x, 1+size(bl,2)+(1:length(m))];
                all_y = [all_y, m];
                all_y_std = [all_y_std, s];
                
                % WASHOUT
                m = zeros(1,size(wo,2));
                s = zeros(2,size(wo,2));
                for k = 1:size(wo,2)
                    m(k) = nanmedian(wo(~bad_inds(:,wo_inds(k)),k),1);
                    if do_conf_int
                        temp = wo(~bad_inds(:,wo_inds(k)),k);
                        bs = zeros(1,num_bootstraps);
                        for z = 1:num_bootstraps
                            bs(z) = nanmedian(temp(randi(length(temp),length(temp),1)),1);
                        end
                        s(:,k) = prctile(bs,[2.5,97.5])';
                    else
                        s(:,k) = [m(k)-nanstd(wo(~bad_inds(:,wo_inds(k)),k),1)./sqrt(size(wo(~bad_inds(:,wo_inds(k)),k),1)); ...
                            m(k)+nanstd(wo(~bad_inds(:,wo_inds(k)),k),1)./sqrt(size(wo(~bad_inds(:,wo_inds(k)),k),1))];
                    end
                end
                all_x = [all_x, 1+size(bl,2)+size(ad,2)+(1:length(m))];
                all_y = [all_y, m];
                all_y_std = [all_y_std, s];
                
                all_x = all_x - 1;
                
                min_y = min([min_y, all_y_std(1,:)]);
                max_y = max([max_y, all_y_std(2,:)]);
                
                hold all;
                patch([all_x, fliplr(all_x)],[all_y_std(1,:), fliplr(all_y_std(2,:))],plot_colors(i,:),'FaceAlpha',0.3,'EdgeAlpha',0.5,'EdgeColor',plot_colors(i,:));
%               plot(all_x,all_y,'o','LineWidth',3,'Color',plot_colors(i,:));
                h(i) = plot(all_x,all_y,'-','LineWidth',3,'Color',plot_colors(i,:));
                
                
                
            case 4
                %% HISTOGRAM DIFFERENCE VERSION
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                d_vaf = all_m - repmat(mean(all_m_bl,2),1,size(all_m,2));
                
                if do_norm
                    if use_cv
                        d_vaf = d_vaf ./ repmat(abs(mean(all_m_cv,2)),1,size(all_m,2));
                    else
                        d_vaf = d_vaf ./ repmat(abs(mean(all_m_bl,2)),1,size(all_m,2));
                    end
                end
                
                bad_inds = zeros(size(d_vaf));
                if remove_outliers
                    for k = 1:size(d_vaf,2)
                        temp = prctile(d_vaf(:,k),[outlier_alpha/2, 100-outlier_alpha/2]);
                        bad_inds(:,k) = d_vaf(:,k) < temp(1) | d_vaf(:,k) > temp(2);
                    end
                end
                
                subplot1(i);
                hold all;
                [N,X] = hist(d_vaf(~bad_inds(:,ad_inds(end)),ad_inds(end)),nbins);
                N = N./sum(N);
                bar(X,N,1,'FaceColor','r','FaceAlpha',0.7);
                [N,X] = hist(d_vaf(~bad_inds(:,wo_inds(end)),wo_inds(end)),nbins);
                N = N./sum(N);
                bar(X,N,1,'FaceAlpha',0.7);
                axis('tight'); V=axis;
                plot([0 0],V(3:4),'k--','LineWidth',2);
                set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[xmin xmax]);
                if i==1
                    ylabel('Density','FontSize',14);
                end
                xlabel('Change in Relative Pseudo-R2','FontSize',14);
                title(array_pairs{i},'FontSize',16);
                max_y = max([max_y, V(4)]);
                
                % do some statistical tests
                [~,p] = ttest2(d_vaf(~bad_inds(:,ad_inds(end)),ad_inds(end)),d_vaf(~bad_inds(:,wo_inds(end)),wo_inds(end)),'tail','left');
                disp([pert{idx_pert} ' ' array_pairs{i} ' - AD2 on WO2: ' num2str(p) '; Mean: ' num2str((mean(d_vaf(~bad_inds(:,wo_inds(end)),wo_inds(end)))-mean(d_vaf(~bad_inds(:,ad_inds(end)),ad_inds(end))))/mean(d_vaf(~bad_inds(:,ad_inds(end)),ad_inds(end))))]);
                [~,p] = ttest2(d_vaf(~bad_inds(:,ad_inds(end)),ad_inds(end)),d_vaf(~bad_inds(:,wo_inds(1)),wo_inds(1)),'tail','left');
                disp([pert{idx_pert} ' ' array_pairs{i} ' - AD2 on WO1: ' num2str(p)]);
                [~,p] = ttest2(d_vaf(~bad_inds(:,ad_inds(end)),ad_inds(end)),d_vaf(~bad_inds(:,ad_inds(1)),ad_inds(1)),'tail','left');
                disp([pert{idx_pert} ' ' array_pairs{i} ' - AD1 on AD2: ' num2str(p) '; Mean: ' num2str((mean(d_vaf(:,ad_inds(end)))-mean(d_vaf(:,ad_inds(1))))/mean(d_vaf(:,ad_inds(1))))]);
                [~,p] = ttest2(d_vaf(~bad_inds(:,wo_inds(end)),wo_inds(end)),d_vaf(~bad_inds(:,ad_inds(1)),ad_inds(1)),'tail','left');
                disp([pert{idx_pert} ' ' array_pairs{i} ' - AD1 on WO2: ' num2str(p)]);
                
                [~,p] = ttest2(d_vaf(~bad_inds(:,wo_inds(end)),wo_inds(end)),d_vaf(~bad_inds(:,bl_inds),bl_inds),'tail','left');
                disp([pert{idx_pert} ' ' array_pairs{i} ' - BL on WO2: ' num2str(p)]);
        end
    end
    
    % some end-of-script tweaking of stuff that is option specific
    switch do_plot
        case 2
            for i = 1:length(array_pairs)
                subplot(2,length(array_pairs),i)
                set(gca,'YLim',[0 max_y]);
                plot([0 0],[0 max_y],'k--','LineWidth',2);
                subplot(2,length(array_pairs),i+length(array_pairs))
                set(gca,'YLim',[0 max_y]);
                plot([0 0],[0 max_y],'k--','LineWidth',2);
            end
        case 3
            
            min_y = -1.5; max_y = 0.25;
            
            for i = 1:length(array_pairs)
                set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 length(ad_inds)+length(wo_inds)],'YLim',[min_y, max_y], ...
                    'XTick',[floor(0+length(ad_inds)/2), ceil(1+length(ad_inds)+length(wo_inds)/2)],'XTickLabel',{'Force','Wash'},'FontSize',14);
                plot([length(bl_inds), length(bl_inds)]+0.5,[min_y, max_y],'k--','LineWidth',2);
                plot(0.5+length(bl_inds)+[length(ad_inds), length(ad_inds)],[min_y, max_y],'k--','LineWidth',2);
            end
            if idx_pert == 1
                ylabel('Norm Change in Rel-Pseudo-R^2','FontSize',14);
            end
            if idx_pert == length(pert)
                legend(h,array_pairs,'Box','off','FontSize',14,'Location','NorthEast');
            end
        case 4
            for i = 1:length(array_pairs)
                subplot1(i);
                if i==length(array_pairs)
                    legend({'Force','Wash'},'FontSize',14);
                end
                set(gca,'YLim',[0 max_y]);
                plot([0 0],[0 max_y],'k--','LineWidth',2);
            end
    end
end