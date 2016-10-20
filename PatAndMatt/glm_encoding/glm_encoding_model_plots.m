%%
clear
clc;
% close all;

dataSummary;

outputSubdir = 'plan';
sessions = {'Chewie','2016-09-15'};
% sessions = {'Mihili','2014-02-17'};

% array_pairs = {'M1-M1','PMd-PMd','PMd-M1'};
array_pairs = {'PMd-M1'};

min_pr2 = 0;
which_metric = 'pr2_full';

do_norm = true;
remove_outliers =false; % supported for 5 and 7 right now
outlier_alpha = 5;
do_conf_int = true; % bootstrapped 95% confidence bounds for error bars
num_bootstraps = 1000;

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

session_idx = strcmpi(filedb.Monkey,sessions(:,1)) & strcmpi(filedb.Date,sessions(:,2));
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
        
        [all_m, all_m_bl,all_m_cv] = deal([]);
        % loop along files and group together
        [total_cells,total_significant] = deal(0);
        
        for file = 1:length(filenames) % loop along sessions
            load(fullfile(rootDir,TDDir,outputSubdir,[perts{idx_pert} '-' array_pairs{model} '_' filenames{file} '.mat']),'results','params');
            
            % parse out relevant params
            bl_inds = find(strcmpi(params.test_epochs,'BL'));
            ad_inds = find(strcmpi(params.test_epochs,'AD'));
            wo_inds = find(strcmpi(params.test_epochs,'WO'));
            
%             temp_pr2 = max(results.pr2_full_cv(:,:,1),[],2);
%             temp_rpr2 = max(results.rpr2_cv(:,:,1),[],2);
%             temp_pr2 = mean(mean(results.pr2_full_cv,3),2);
%             temp_rpr2 = mean(mean(results.rpr2_cv,3),2);
%             good_idx = temp_pr2(:,1) > min_pr2_full & temp_rpr2(:,1) > min_rpr2;

            temp_pr2 = mean(mean(results.([which_metric '_cv']),3),2);
            good_idx = temp_pr2(:,1) > min_pr2;

            temp_metric = results.(which_metric);
            temp_metric = temp_metric(good_idx,:,:);
            
            % get some stats
            total_significant = total_significant+sum(good_idx);
            total_cells = total_cells+length(good_idx);
            
            if ~isempty(bl_inds)
                temp_metric_bl = squeeze(mean(results.(which_metric)(:,bl_inds,:),2));
                temp_metric_bl = temp_metric_bl(good_idx,:,:);
            end
            
            temp_metric_cv = results.([which_metric '_cv']);
            
            % if this file is
            
            all_m_cv = [all_m_cv; mean(temp_metric_cv(good_idx,:,:),3)];
            % pool together metrics from all files for plotting
            all_m = [all_m; mean(temp_metric,3)];
        end
        
        disp([array_pairs{model} ' - % cells with significant rel-pseudo-r2: ' num2str(total_significant) '/' num2str(total_cells)]);
        
        disp([array_pairs{model} ' - Mean baseline metric: ' num2str(mean(mean(all_m_cv,2))) ' +/- ' num2str(std(mean(all_m_cv,2)))]);
        
        all_shit = [all_shit; mean(all_m_bl,2)];
        
        % OVER TIME
        plot_inds = [bl_inds;ad_inds;wo_inds];
        
        d_vaf = all_m - repmat(mean(all_m_cv,2),1,size(all_m,2));
        
        if do_norm
            d_vaf = d_vaf ./ repmat(abs(mean(all_m_cv,2)),1,size(all_m,2));
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
        
%         min_y = min([min_y, all_y_std(1,:)]);
%         max_y = max([max_y, all_y_std(2,:)]);
        
        hold all;
        patch([all_x, fliplr(all_x)],[all_y_std(1,:), fliplr(all_y_std(2,:))],plot_colors(model,:),'FaceAlpha',0.3,'EdgeAlpha',0.5,'EdgeColor',plot_colors(model,:));
        h(model) = plot(all_x,all_y,'+','LineWidth',3,'Color',plot_colors(model,:));
%         h(model) = plot(all_x,all_y,'-','LineWidth',3,'Color',plot_colors(model,:));
        
    end
    
    % some end-of-script tweaking of stuff    
    for model = 1:length(array_pairs)
        set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 length(ad_inds)+length(wo_inds)],'YLim',[min_y, max_y], ...
            'XTick',[floor(0+length(ad_inds)/2), ceil(1+length(ad_inds)+length(wo_inds)/2)],'XTickLabel',{'Force','Wash'},'FontSize',14);
        plot([length(bl_inds), length(bl_inds)]+0.5,[min_y, max_y],'k--','LineWidth',2);
        plot(0.5+length(bl_inds)+[length(ad_inds), length(ad_inds)],[min_y, max_y],'k--','LineWidth',2);
    end
    if idx_pert == 1
        ylabel('Norm Change in Rel-Pseudo-R^2','FontSize',14);
    end
    if idx_pert == length(perts)
        legend(h,array_pairs,'Box','off','FontSize',14,'Location','NorthEast');
    end
    
end