clear;
clc;
close all;

dataSummary;


params_file = '';% optionally give .mat file containing params to use
% Must be in the outputSubdir

sessions = { ...
    'Chewie','2016-09-09'; ... % VR
    'Chewie','2016-09-12'; ...
    'Chewie','2016-09-14'; ...
    'Chewie','2016-10-06'; ...
%     'Mihili','2014-03-03'; ...
%     'Mihili','2014-03-04'; ...
%     'Mihili','2014-03-06'; ...
    'Chewie','2016-09-15'; ... % CF
%     'Chewie','2016-09-19'; ...
    'Chewie','2016-10-05'; ...
    'Chewie','2016-10-07'; ...
    'Chewie','2016-10-11'; ...
%         'Mihili','2014-02-03'; ...
%         'Mihili','2014-02-17'; ...
%         'Mihili','2014-02-18'; ...
%         'Mihili','2014-03-07'; ...
    };


monkeys = unique(sessions(:,1));
tasks = {'CO'};
pert = 'FF';
dates = sessions(:,2);

ref_epoch = {'BL',[0 1]};
test_epochs = {'AD',[0 0.5]; ...
    'WO',[0 0.5]};

mahal_arrays = {'M1','PMd','potent','null'};
use_cues = {'prego','go'};%tgt,prego,go,move,peak

%%
if isempty(params_file)
    dt                   = 0.01; % time step size for data
    bin_size             = 5;    % how many samples to group together when rebinning (bin width is num_samples*dt)
    
    result_codes         = {'R','I'}; % which trials to include
    
    % how to truncate trials {idx name, number of bins after}
    start_idx      = {'idx_target_on',0}; %{'idx_target_on',-4}
    end_idx        = {'idx_trial_end',0}; %{'idx_go_cue',4}
    %   NOTE: this is after rebinning at the moment
    
    block_size_fr_test   = 100;   % how many trials to group for FR test
    fr_test_alpha        = 1e-3; % p value cut off for t-test
    fr_min               = 0.15; % minimum session-wide spiking for inclusion
    
    pca_dims             = struct('M1',1:8,'PMd',1:16); % ...'ARRAY','all' or specify which dimensions, e.g. 1:30...
    
else % load up params from whatever file
    load(params_file,'cv_params');
    fn = fieldnames(params);
    for i = 1:length(fn)
        eval([fn{i} ' = cvparams.' fn{i} ';']);
    end, clear i good_cells cov_array pred_array pert epochs pca_w;
end


plot_colors = [0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840; ...
    0 0 1];
%% Here is where I check all of the inputs

%%
disp(['Starting perturbation ' pert '...']);

use_date_idx = ismember(filedb.Monkey,monkeys) & ismember(filedb.Perturbation,pert) & ismember(filedb.Task,tasks);
if ~isempty(dates)
    use_date_idx = use_date_idx & ismember(filedb.Date,dates);
end
use_files = find(use_date_idx);

[all_dist_ad,all_dist_ad] = deal(cell(length(mahal_arrays),length(use_files),length(use_cues)));

h = figure;
for idx_file = 1:length(use_files)
    disp(['File ' num2str(idx_file) ' of ' num2str(length(use_files)) '...'])
    
    epochs = filedb.Epochs{use_files(idx_file)};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Package up parameters for saving
    clear params;
    params.block_size_fr_test = block_size_fr_test;
    params.fr_test_alpha = fr_test_alpha;
    params.fr_min = fr_min;
    params.dt = dt;
    params.bin_size = bin_size;
    params.pca_dims = pca_dims;
    params.result_codes = result_codes;
    params.start_idx = start_idx;
    params.end_idx = end_idx;
    cv_params = params; % these will be saved separately
    params.pert = pert;
    params.filedb = filedb;
    params.epochs = epochs;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load the trial data file
    filename = [filedb.Monkey{use_files(idx_file)} '_' filedb.Task{use_files(idx_file)} '_' filedb.Perturbation{use_files(idx_file)} '_' filedb.Date{use_files(idx_file)}];
    load(fullfile(rootDir,TDDir,[filename '.mat']));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % filter out neurons
    trial_data = getCommonUnits(trial_data);
    
    arrays = {'M1','PMd'};
    good_cells = cell(1,length(arrays));
    for array = 1:length(arrays);
        % make sure firing is significantly non-zero across the
        % whole session. Usually this means there is a sort problem
        % or there is noise and I lose a neuron, a real neuron
        % wouldn't completely shut off.
        all_fr=cell2mat(cellfun(@(x) sqrt(sum(x,1))', ...
            cellfun(@(x) full(x),{trial_data.([arrays{array} '_spikes'])}, ...
            'Uni',0),'Uni',0));
        
        % additionally, make sure the baseline firing (with the
        % model fit) is sufficiently high. I do this separately
        % from the last one because I want to allow neurons to
        % increase or decrease firing as they see fit during
        % learning
        bl_fr = [];
        trial_idx = find(strcmpi({trial_data.epoch},'BL'));
        for iTrial = trial_idx
            idx = trial_data(iTrial).(start_idx{1})+start_idx{2}:trial_data(iTrial).(end_idx{1})+end_idx{2};
            temp = trial_data(iTrial).([arrays{array} '_spikes']);
            bl_fr = [bl_fr, (sqrt(sum(temp(idx,:),1))./(size(temp,1)*dt))'];
        end
        
        % ensure all blocks are significantly non-zero
        blocks = 1:block_size_fr_test:size(all_fr,2);
        p = zeros(size(all_fr,1),length(blocks));
        for block = 2:length(blocks)
            for unit = 1:size(all_fr,1)
                [~,p(unit,block)] = ttest(all_fr(unit,blocks(block-1):blocks(block)),0,'tail','right');
            end
        end
        good_cells{array} = find(all(p <= fr_test_alpha,2) & mean(bl_fr,2) >= fr_min)';
        
        disp(['Removing ' num2str(size(all_fr,1) - length(good_cells{array})) ' low-firing cells...']);
        
        % take bad cells out of trial_data
        for trial = 1:length(trial_data)
            temp = trial_data(trial).([arrays{array} '_spikes']);
            trial_data(trial).([arrays{array} '_spikes']) = temp(:,good_cells{array});
            temp = trial_data(trial).([arrays{array} '_unit_guide']);
            trial_data(trial).([arrays{array} '_unit_guide']) = temp(good_cells{array},:);
        end
    end, clear all_fr bl_fr num_blocks i array r s temp blocks;
    params.good_cells = good_cells;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Re-bin data
    if bin_size > 1
        disp('Binning...');
        trial_data = truncateAndBin(trial_data,bin_size);
        params.dt = dt*bin_size;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute PCA
    disp('Getting potent/null space.')
    % do SVD to get null/potent spaces
    [V_potent, V_null, w_pmd, w_m1] = getPotentSpace(trial_data,'PMd','M1',params);
    params.w_m1 = w_m1;
    params.w_pmd = w_pmd;
    params.V_potent = V_potent;
    params.V_null = V_null;
    
    
    
    %%
    % loop along training trials and average together
    
    utheta = unique([trial_data.target_direction]);
    
    w.M1 = w_m1; w.PMd = w_pmd;
    
    for c = 1:length(use_cues)
        switch lower(use_cues{c})
            case 'tgt'
                which_cue = 'target_on';
                idx_lag = 2;
                bins_before = 0;
                bins_after = 4;
            case 'prego'
                which_cue = 'go_cue';
                idx_lag = 0;
                bins_before = 4;
                bins_after = 0;
            case 'go'
                which_cue = 'go_cue';
                idx_lag = 0;
                bins_before = 0;
                bins_after = 4;
            case 'move'
                which_cue = 'movement_on';
                idx_lag = 0;
                bins_before = 1;
                bins_after = 3;
            case 'peak'
                which_cue = 'peak_speed';
                idx_lag = 0;
                bins_before = 2;
                bins_after = 2;
        end
        
        for a = 1:length(mahal_arrays)
            %subplot(1,length(arrays),a); hold all;
            array = mahal_arrays{a};
            
            num_dims = max(pca_dims.M1);
            
            bl_idx = find(getTDidx(trial_data,'epoch',ref_epoch{1}));
            bl_idx = bl_idx(1+floor(ref_epoch{2}(1)*length(bl_idx)):floor(ref_epoch{2}(2)*length(bl_idx)));
            ad_idx = find(getTDidx(trial_data,'epoch',test_epochs{1,1}));
            ad_idx = ad_idx(1+floor(test_epochs{1,2}(1)*length(ad_idx)):floor(test_epochs{1,2}(2)*length(ad_idx)));
            wo_idx = find(getTDidx(trial_data,'epoch',test_epochs{2,1}));
            wo_idx = wo_idx(1+floor(test_epochs{2,2}(1)*length(wo_idx)):floor(test_epochs{2,2}(2)*length(wo_idx)));
            
            bl_proj = cell(1,length(utheta));
            for th = 1:length(utheta)
                
                trial_idx = find([trial_data(bl_idx).target_direction]==utheta(th));
                
                bl = zeros(1+bins_before+bins_after,num_dims,length(trial_idx));
                for i = 1:length(trial_idx)
                    trial = trial_idx(i);
                    
                    switch lower(array)
                        case 'potent'
                            temp = sqrt(trial_data(trial).PMd_spikes);
                            w.PMd;
                            temp = [ones(size(temp,1),1) temp(:,pca_dims.PMd)]*V_potent;
                        case 'null'
                            temp = sqrt(trial_data(trial).PMd_spikes)*w.PMd;
                            temp = [ones(size(temp,1),1) temp(:,pca_dims.PMd)]*V_null;
                        otherwise
                            temp = sqrt(trial_data(trial).([array '_spikes']))*w.(array);
                    end
                    idx = trial_data(trial).(['idx_' which_cue])-bins_before:trial_data(trial).(['idx_' which_cue])+bins_after;
                    idx = idx + idx_lag;
                    bl(:,:,i) = temp(idx,1:num_dims);
                end
                bl_proj{th} = squeeze(mean(bl,1));
            end
            
            % loop along ad trials and get distance at idx from bl
            [mahal_dist_ad, mahal_dist_wo] = deal(zeros(1,length(ad_idx)));
            for i = 1:length(ad_idx)
                trial = ad_idx(i);
                bl = bl_proj{utheta == trial_data(trial).target_direction};
                
                switch lower(array)
                    case 'potent'
                        temp = sqrt(trial_data(trial).PMd_spikes)*w.PMd;
                        temp = [ones(size(temp,1),1) temp(:,pca_dims.PMd)]*V_potent;
                    case 'null'
                        temp = sqrt(trial_data(trial).PMd_spikes)*w.PMd;
                        temp = [ones(size(temp,1),1) temp(:,pca_dims.PMd)]*V_null;
                    otherwise
                        temp = sqrt(trial_data(trial).([array '_spikes']))*w.(array);
                end
                
                idx = trial_data(trial).(['idx_' which_cue])-bins_before:trial_data(trial).(['idx_' which_cue])+bins_after;
                idx = idx +idx_lag;
                ad_proj = squeeze(mean(temp(idx,1:num_dims),1));
                
                mahal_dist_ad(i) = mahal(ad_proj,bl');
            end
            
            % loop along wo trials and get distance at idx from bl
            for i = 1:length(wo_idx)
                trial = wo_idx(i);
                bl = bl_proj{utheta == trial_data(trial).target_direction};
                
                switch lower(array)
                    case 'potent'
                        temp = sqrt(trial_data(trial).PMd_spikes)*w.PMd;
                        temp = [ones(size(temp,1),1) temp(:,pca_dims.PMd)]*V_potent;
                    case 'null'
                        temp = sqrt(trial_data(trial).PMd_spikes)*w.PMd;
                        temp = [ones(size(temp,1),1) temp(:,pca_dims.PMd)]*V_null;
                    otherwise
                        temp = sqrt(trial_data(trial).([array '_spikes']))*w.(array);
                end
                
                idx = trial_data(trial).(['idx_' which_cue])-bins_before:trial_data(trial).(['idx_' which_cue])+bins_after;
                idx = idx +idx_lag;
                ad_proj = squeeze(mean(temp(idx,1:num_dims),1));
                
                mahal_dist_wo(i) = mahal(ad_proj,bl');
            end
            
            all_dist_ad{a,idx_file,c} = mahal_dist_ad;
            all_dist_wo{a,idx_file,c} = mahal_dist_wo;
        end
    end
    
    
end % end file loop


%% Plot some shit

bin_size = 1;

close all;
figure('Position',[300 25 900 950]);
subplot1(length(use_cues),length(mahal_arrays),'Gap',[0.02, 0.01]);
save_p = zeros(1,length(use_cues)*length(mahal_arrays));

    ymin = Inf;
    ymax = -Inf;
plot_order = [1,2,5,6,3,4,7,8];
for c = 1:length(use_cues)

    for a = 1:length(mahal_arrays)
        h((c-1)*length(mahal_arrays) + a) = subplot1(plot_order((c-1)*length(mahal_arrays) + a)); hold all;
        
        d = all_dist_ad(a,:,c);
        
        d = cell2mat(cellfun(@(x) sqrt(x(1:min(cellfun(@(x) length(x),d)))),d,'Uni',0)');
        
        for i = 1:size(d,1)
            [~,~,~,rint,~] = regress(d(i,:)',[ones(size(d,2),1), (1:size(d,2))']);
            d(i,rint(:,1)>0) = NaN;
        end
        
        d = d ./ repmat(nanmean(d,2),1,size(d,2));
        
        bins = 1:bin_size:size(d,2);
        [save_m,save_s] = deal(zeros(1,length(bins)-1));
        save_all = cell(1,length(bins)-1);
        for i = 1:length(bins)-1
            temp = d(:,1+bins(i):bins(i+1));
            if 1 % pool
                temp = reshape(temp,numel(temp),1);
                temp(isnan(temp)) = [];
            else
                temp = nanmean(temp,2);
            end
            
            save_m(i) = median(temp,1);
            save_s(i) = std(temp,[],1)/sqrt(length(temp));
            save_all{i} = temp';
        end
        
        [~,~,~,rint,~] = regress(save_m',[ones(length(save_m),1), (1:length(save_m))']);
        
        % identify outliers from linear fit
        bad_idx = rint(:,1) > 0;
        save_m(bad_idx) = [];
        save_s(bad_idx) = [];
        save_all(bad_idx) = [];
        
        %save_m = save_m - save_m(1);
        
        [b,bint,r,rint,s] = regress(save_m',[ones(length(save_m),1), (1:length(save_m))']);
        
        save_p((c-1)*length(mahal_arrays) + a) = s(3);
        
        ymin = min([ymin,min(save_m-save_s)]);
        ymax = max([ymin,max(save_m+save_s)]);
        
        plot((1:length(save_m))*bin_size,save_m,'ko');
        plot([1:length(save_m); 1:length(save_m)]*bin_size,[save_m-save_s; save_m+save_s],'k-');
        
        if c == 1
            plot((1:length(save_m))*bin_size,b(1)+b(2)*(1:length(save_m)),'b-','LineWidth',2);
        else 
            plot((1:length(save_m))*bin_size,b(1)+b(2)*(1:length(save_m)),'r-','LineWidth',2);
        end
        
        axis('tight');
        set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 length(bins)]*bin_size);
        
        if plot_order((c-1)*length(mahal_arrays) + a) == 1 || plot_order((c-1)*length(mahal_arrays) + a) == length(mahal_arrays)+1
            ylabel('Norm Change in Mahal Dist','FontSize',14);
        end
        if plot_order((c-1)*length(mahal_arrays) + a) > length(mahal_arrays)
            xlabel('Trials','FontSize',14);
        end
        if a==1
            if c==1
                title('Planning Period','FontSize',16,'Color','b');
            elseif c==2
                title('Reaction Period','FontSize',16,'Color','r');
            end
        end
        
        g = [];
        for i = 1:length(save_all)
            g = [g, i*ones(1,length(save_all{i}))];
        end
        p = anovan(cell2mat(save_all)',g','display','off');
        disp([use_cues{c} ' - ' mahal_arrays{a} ' - p = ' num2str(p,4)]);
    end
    
    disp(' ');
end
for c = 1:length(use_cues)
    for a = 1:length(mahal_arrays)
        set(h((c-1)*length(mahal_arrays) + a),'YLim',[ymin,ymax]);
        subplot1(plot_order((c-1)*length(mahal_arrays) + a));
        text(0.2*bin_size,ymax-0.05*(ymax-ymin),[mahal_arrays{a} ' p=' num2str(save_p((c-1)*length(mahal_arrays) + a),3)],'FontSize',14);
    end
end
linkaxes(h);


