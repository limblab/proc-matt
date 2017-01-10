% plot trajectories in some space
clear;
clc;
close all;

dataSummary;

monkey = 'Chewie';
date = '2016-10-05';

use_file = find(ismember(filedb.Monkey,monkey) & ismember(filedb.Date,date));


plot_colors = [0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840; ...
    0 0 1];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the trial data file
filename = [filedb.Monkey{use_file} '_' filedb.Task{use_file} '_' filedb.Perturbation{use_file} '_' filedb.Date{use_file}];
load(fullfile(rootDir,TDDir,[filename '.mat']));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params_file = ['F:\TrialDataFiles\trainad_potent_all\FF-PMd-M1_Chewie_CO_FF_2016-09-15_cv.mat'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the trial data file
load(params_file,'cv_params');
fn = fieldnames(cv_params);
for i = 1:length(fn)
    eval([fn{i} ' = cv_params.' fn{i} ';']);
end, clear i good_cells cov_array pred_array pert epochs pca_w;
params = cv_params;
params.arrays = {'M1','PMd'};
params.bin_size = 2;

params.do_rcb = false;
params.kin_signals = {};
params.fr_min = 0.3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[trial_data, params] = glm_process_trial_data(trial_data,params);

trial_data = trial_data(ismember({trial_data.result},{'R'}));

%%
[V_potent, V_null, w_pmd, w_m1] = getPotentSpace(trial_data,'PMd','M1',params);

% add projections to trial data
for trial = 1:length(trial_data)
    temp = sqrt(trial_data(trial).M1_spikes);
    temp = smoothSpikesForPCA(temp,0.02,0.04);
    trial_data(trial).M1_pca = temp*w_m1;
    
    temp = sqrt(trial_data(trial).PMd_spikes);
    temp = smoothSpikesForPCA(temp,0.02,0.04);
    trial_data(trial).PMd_pca = temp*w_pmd;
    
    temp = temp*w_pmd;
    trial_data(trial).Potent_pca = temp(:,pca_dims.PMd)*V_potent;
    trial_data(trial).Null_pca = temp(:,pca_dims.PMd)*V_null;
end

%%
% get GPFA trajectories
% try the three different variants to make sure they work
params.bin_w = 20;
params.data_bin_w = 20;
params.xdim = 6;
params.kernsd = 30;
params.arrays = 'M1';
[trial_data, M1_gpfa] = run_gpfa(trial_data,params);

%%
params.xdim = 12;
params.arrays = 'PMd';
[trial_data, PMd_gpfa] = run_gpfa(trial_data,params);

%%
% % % for trial = 1:length(trial_data)
% % %     trial_data(trial).M1_gpfa = M1_gpfa.trajectories(trial).xsm';
% % %     trial_data(trial).PMd_gpfa = PMd_gpfa.trajectories(trial).xsm';
% % % end
%% Find GPFA potent/null spaces
% get GPFA matrix
m1_dims = 1:4;
pmd_dims = 1:8;
[x,y] = deal([]);
for trial = 1:length(trial_data)
    y = [y; trial_data(trial).M1_gpfa(:,m1_dims)];
    x = [x; trial_data(trial).PMd_gpfa(:,pmd_dims)];
end

% find the model
W = zeros( size(y,2), size(x,2) );
for i = 1:size(y,2)
    [b_pc, ~, ~, ~, stats_this] = regress(y(:,i),x);
    % fill MIMO matrix W
    W(i,:) = b_pc';
end
% do SVD of weights
[U, S, V]                   = svd( W );
% The output potent spaces is defined by the first m columns of V', where m
% is the number of dimensions of the output
V_potent                    = V(1:size(y,2),:)';
V_null                      = V(size(y,2)+1:end,:)';

for trial = 1:length(trial_data)
    trial_data(trial).Potent_gpfa = trial_data(trial).PMd_gpfa(:,pmd_dims)*V_potent;
    trial_data(trial).Null_gpfa = trial_data(trial).PMd_gpfa(:,pmd_dims)*V_null;
end


%% plot mahal distance heat map
clc
close all;

use_aligns = {'idx_go_cue',[10,0], 'Planning'; ...
    'idx_go_cue',[-10,20], 'Reaction'};
use_dims = 1:8;
bin_size = 20;

% arrays = {'M1_gpfa','PMd_gpfa'};
arrays = {'M1_pca','PMd_pca'};

utheta = unique([trial_data.target_direction]);


[mahal_dist_ad,mahal_dist_wo,euclid_dist_ad,euclid_dist_wo,cue_dist] = deal(cell(length(arrays),size(use_aligns,1)));
for a = 1:length(arrays);
    array = arrays{a};
    for c = 1:size(use_aligns,1)
        align_idx = use_aligns{c,1};
        align_window = use_aligns{c,2};
        
        % get baseline for each direction
        all_bl = cell(1,length(utheta));
        for th = 1:length(utheta)
            bl_idx = find(getTDidx(trial_data,'epoch','bl'));
            
            trial_idx = find([trial_data(bl_idx).target_direction]==utheta(th));
            
            bl = zeros(1+sum(align_window),length(use_dims),length(trial_idx));
            for i = 1:length(trial_idx)
                trial = trial_idx(i);
                temp = trial_data(trial).(array);
                idx = trial_data(trial).(align_idx)-align_window(1):trial_data(trial).(align_idx)+align_window(2);
                bl(:,:,i) = temp(idx,use_dims);
            end
            all_bl{th} = bl;
        end
        
        trial_idx = find(getTDidx(trial_data,'epoch','ad'));
%         trial_idx = trial_idx(1:floor(0.5*length(trial_idx)));
        [mahal_dist,euclid_dist] = deal(zeros(1+sum(align_window),length(trial_idx)));
        for i = 1:length(trial_idx)
            trial = trial_idx(i);
            bl = all_bl{utheta==trial_data(trial).target_direction};
            
            temp = trial_data(trial).(array);
            idx = trial_data(trial).(align_idx)-align_window(1):trial_data(trial).(align_idx)+align_window(2);
            ad = temp(idx,use_dims);
            for j = 1:size(ad,1)
                mahal_dist(j,i) = mahal(ad(j,:),squeeze(bl(j,:,:))');
                euclid_dist(j,i) = norm(ad(j,:)' - mean(squeeze(bl(j,:,:)),2));
            end
        end
        mahal_dist_ad{a,c} = mahal_dist;
        euclid_dist_ad{a,c} = euclid_dist;
        %cue_dist{a,c} = euclid_dist(1+align_window(1),:);
        
        trial_idx = find(getTDidx(trial_data,'epoch','wo'));
%         trial_idx = trial_idx(1:floor(0.5*length(trial_idx)));
        [mahal_dist,euclid_dist] = deal(zeros(1+sum(align_window),length(trial_idx)));
        for i = 1:length(trial_idx)
            trial = trial_idx(i);
            bl = all_bl{utheta==trial_data(trial).target_direction};
            
            temp = trial_data(trial).(array);
            idx = trial_data(trial).(align_idx)-align_window(1):trial_data(trial).(align_idx)+align_window(2);
            ad = temp(idx,use_dims);
            for j = 1:size(ad,1)
                mahal_dist(j,i) = mahal(ad(j,:),squeeze(bl(j,:,:))');
                euclid_dist(j,i) = norm(ad(j,:)' - mean(squeeze(bl(j,:,:)),2));
            end
        end
        mahal_dist_wo{a,c} = mahal_dist;
        euclid_dist_wo{a,c} = euclid_dist;
    end
end


figure;
subplot1(length(arrays),size(use_aligns,1));
for a = 1:length(arrays)
    for c = 1:size(use_aligns,1)
        align_idx = use_aligns{c,1};
        align_window = use_aligns{c,2};
        align_title = use_aligns{c,3};
        
        subplot1((a-1)*size(use_aligns,1)+c); hold all;
        temp = sqrt(mahal_dist_ad{a,c});
        
        % group some trials together
        bins = 1:bin_size:size(temp,2);
        new_temp = zeros(size(temp,1),length(bins)-1);
        for bin = 1:length(bins)-1
            new_temp(:,bin) = mean(temp(:,bins(bin):bins(bin+1)-1),2);
        end
        
        imagesc(new_temp');
        axis('tight');
        V = axis;
        plot([align_window(1),align_window(1)],V(3:4),'w-','LineWidth',2);
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[1,sum(align_window)]);
        
        if a == 1, title(align_title,'FontSize',14); end
        if a == length(arrays),
%             set(gca,'XTick',[1, 1+align_window(1), 1+sum(align_window)],'XTickLabel',{num2str(-align_window(1)),strrep(align_idx(5:end),'_','-'),num2str(align_window(2))});
            xlabel('Time in trial','FontSize',14);
        end
        if c == 1, ylabel(arrays{a},'FontSize',16); end
        
    end
end

% % plot correlation with learning
% for a = 1:length(arrays)
%     figure; subplot1(1,size(use_aligns,1))
%     ymin = Inf; ymax = -Inf;
%     title(arrays{a})
%     for c = 1:size(use_aligns,1)
%         subplot1(c); hold all;
%         temp = euclid_dist_ad{a,c};
%         temp = mean(temp,1);
%         
%         % group some trials together
%         bins = 1:bin_size:size(temp,2);
%         new_temp = zeros(size(temp,1),length(bins)-1);
%         for bin = 1:length(bins)-1
%             new_temp(:,bin) = mean(temp(:,bins(bin):bins(bin+1)-1),2);
%         end
%         plot(new_temp,'+');
%         ymin = min([ymin,min(new_temp)]);
%         ymax = max([ymax,max(new_temp)]);
%         [b,~,~,~,s] = regress(new_temp',[ones(length(new_temp),1), (1:length(new_temp))']);
%         plot(1:length(new_temp),b(1)+b(2)*(1:length(new_temp)),'k-','LineWidth',2);
%         title(['p = ' num2str(s(3),3)]);
%     end
%     subplot1(1); ylabel('Distance in Neural Space','FontSize',14);
%     for c = 1:size(use_aligns,1)
%         subplot1(c);
%         axis('tight');
%         set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[ymin,ymax]);
%         xlabel('Curl Field Trials','FontSize',14);
%     end
% end





%% plot trajectories
close all;

align_idx = {'idx_go_cue','idx_go_cue'};
align_window = [15,10]; % bins before, after
idx_lag = 0;
use_dims = [1,2,3];

arrays = {'PMd_gpfa'};

utheta = unique([trial_data.target_direction]);

for th = 1:length(utheta)
    figure;
    for a = 1:length(arrays);
        subplot(1,length(arrays),a); hold all;
        array = arrays{a};
        
        bl_idx = find(getTDidx(trial_data,'epoch','bl'));
        
        trial_idx = find([trial_data(bl_idx).target_direction]==utheta(th));
        
        bl = zeros(1+align_window(1)+align_window(2),length(use_dims),length(trial_idx));
        for i = 1:length(trial_idx)
            trial = trial_idx(i);
            temp = trial_data(trial).(array);
            idx = trial_data(trial).(align_idx{1})-align_window(1):trial_data(trial).(align_idx{2})+align_window(2);
            idx = idx + idx_lag;
            bl = temp(idx,use_dims);
            
            plot3(bl(:,1),bl(:,2),bl(:,3),'-','LineWidth',1,'Color','k');%plot_colors(th,:));
            plot3(bl(1+align_window(1),1),bl(1+align_window(1),2),bl(1+align_window(1),3),'o','LineWidth',3,'Color','k');%plot_colors(th,:));
        end
        
        bl_idx = find(getTDidx(trial_data,'epoch','ad'));
        trial_idx = find([trial_data(bl_idx).target_direction]==utheta(th));
        trial_idx = trial_idx(1:5);
        
        bl = zeros(1+align_window(1)+align_window(2),length(use_dims),length(trial_idx));
        for i = 1:length(trial_idx)
            trial = trial_idx(i);
            temp = trial_data(trial).(array);
            idx = trial_data(trial).(align_idx{1})-align_window(1):trial_data(trial).(align_idx{2})+align_window(2);
            idx = idx + idx_lag;
            bl = temp(idx,use_dims);
            
            plot3(bl(:,1),bl(:,2),bl(:,3),'-','LineWidth',1,'Color','r');%plot_colors(th,:));
            plot3(bl(1+align_window(1),1),bl(1+align_window(1),2),bl(1+align_window(1),3),'o','LineWidth',3,'Color','r');%plot_colors(th,:));
        end
        
                bl_idx = find(getTDidx(trial_data,'epoch','ad'));
        trial_idx = find([trial_data(bl_idx).target_direction]==utheta(th));
        trial_idx = trial_idx(end-4:end);
        
        bl = zeros(1+align_window(1)+align_window(2),length(use_dims),length(trial_idx));
        for i = 1:length(trial_idx)
            trial = trial_idx(i);
            temp = trial_data(trial).(array);
            idx = trial_data(trial).(align_idx{1})-align_window(1):trial_data(trial).(align_idx{2})+align_window(2);
            idx = idx + idx_lag;
            bl = temp(idx,use_dims);
            
            plot3(bl(:,1),bl(:,2),bl(:,3),'-','LineWidth',1,'Color','b');%plot_colors(th,:));
            plot3(bl(1+align_window(1),1),bl(1+align_window(1),2),bl(1+align_window(1),3),'o','LineWidth',3,'Color','b');%plot_colors(th,:));
        end
        
        
        title(strrep(array,'_','-'),'FontSize',16);
        axis('square');
        
    end
end



