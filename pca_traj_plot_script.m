% load data
clear;
close all;
clc;

filenames = { ...
    '/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-07.mat', ...
    '/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-11.mat', ...
    '/Users/mattperich/Data/TrialDataFiles/Chewie_CO_FF_2016-10-13.mat', ...
    };

trial_data = loadTDfiles(filenames,{'truncateAndBin',1},{'removeBadNeurons',struct('min_fr',1)});
trial_data = pruneBadTrials(td,struct('ranges', ...
    {{'idx_go_cue','idx_movement_on',[5,50]}}));

%%
trial_data = smoothSpikes(trial_data,struct('do_smoothing',true,'sqrt_transform',true,'kernel_SD',0.03));
align_idx = 'idx_target_on';
t_before = 10;
t_after = 30;
td1 = truncateAndBin(trial_data,{align_idx,-t_before},{align_idx,t_after});
align_idx = 'idx_go_cue';
t_before = 20;
t_after = 20;
td2 = truncateAndBin(trial_data,{align_idx,-t_before},{align_idx,t_after});
align_idx = 'idx_movement_on';
t_before = 20;
t_after = 50;
td3 = truncateAndBin(trial_data,{align_idx,-t_before},{align_idx,t_after});
td = appendTDs(td1,td2,td3);

%%
% get subspaces
num_dims = 8;
pca_params = struct( ...
    'in_array','PMd', ...
    'out_array','M1', ...
    'in_dims',2*num_dims, ...
    'out_dims',num_dims, ...
    'sqrt_transform',false,...
    'do_smoothing', false, ...
    'bin_size', 0.01, ...
    'kernel_SD',0.02, ...
    'trial_avg',true, ...
    'trial_avg_cond',{{'target_direction','epoch','date'}},...
    'do_plot',false);

dates = unique({td.date});

new_td = [];
for i = 1:length(dates)
    [~,temp_td] = getTDidx(td,'date',dates{i});
    [temp_td,temp] = getPotentSpace(temp_td,pca_params);
    temp_td = getPCA(temp_td,temp.w_out,temp.mu_out,struct('array','M1','do_smoothing',false));
    new_td = [new_td, temp_td];
end
td = new_td;

%%
close all;
clc;
num_outs = 10;
do_corr = false;

dates = unique({td.date});

[all_c_bl,all_c_ad] = deal(zeros(num_outs,2));
[all_c_bl_file,all_c_ad_file] = deal(zeros(num_outs,2,length(dates)));
figure;
for i = 1:num_outs
    which_one = 'potent';
    
    [pred1, pred2, out1, out2] = deal([]);
    for file = 1:length(dates)
        [~,temp_td] = getTDidx(td,'epoch','BL','date',dates{file});
        in_var = cat(1,temp_td.(which_one));
        in_var = [ones(size(in_var,1),1), in_var(:,1:num_dims)];
        out_var = cat(1,temp_td.M1_pca);
        out_var = out_var(:,i);
        
        b = in_var\out_var;
        pred1 = [pred1; in_var*b];
        out1 = [out1; out_var];
        
        if do_corr
            c = corrcoef(out_var,in_var*b); c = c(1,2);
        else
            c = CalculateR2(out_var,in_var*b);
        end     
        all_c_bl_file(i,1,file) = c;
        
        [~,temp_td] = getTDidx(td,'epoch','AD','date',dates{file});
        in_var = cat(1,temp_td.(which_one));
        in_var = [ones(size(in_var,1),1), in_var(:,1:num_dims)];
        out_var = cat(1,temp_td.M1_pca);
        out_var = out_var(:,i);
        pred2 = [pred2; in_var*b];
        out2 = [out2; out_var];
        
        if do_corr
            c = corrcoef(out_var,in_var*b); c = c(1,2);
        else
            c = CalculateR2(out_var,in_var*b);
        end     
        all_c_ad_file(i,1,file) = c;
    end
    pred_bl_p = pred1;
    pred_ad_p = pred2;
    out_bl = out1;
    out_ad = out2;
    
    %%%%%%%%%%%%%%%%%
    which_one = 'null';
    
    [pred1, pred2] = deal([]);
    for file = 1:length(dates)
        [~,temp_td] = getTDidx(td,'epoch','BL','date',dates{file});
        in_var = cat(1,temp_td.(which_one));
        in_var = [ones(size(in_var,1),1), in_var(:,1:num_dims)];
        out_var = cat(1,temp_td.M1_pca);
        out_var = out_var(:,i);
        
        b = in_var\out_var;
        pred1 = [pred1; in_var*b];
        
        if do_corr
            c = corrcoef(out_var,in_var*b); c = c(1,2);
        else
            c = CalculateR2(out_var,in_var*b);
        end     
        all_c_bl_file(i,2,file) = c;
        
        [~,temp_td] = getTDidx(td,'epoch','AD','date',dates{file});
        in_var = cat(1,temp_td.(which_one));
        in_var = [ones(size(in_var,1),1), in_var(:,1:num_dims)];
        pred2 = [pred2; in_var*b];
        out_var = cat(1,temp_td.M1_pca);
        out_var = out_var(:,i);
        
        if do_corr
            c = corrcoef(out_var,in_var*b); c = c(1,2);
        else
            c = CalculateR2(out_var,in_var*b);
        end        
        all_c_ad_file(i,2,file) = c;
    end
    pred_bl_n = pred1;
    pred_ad_n = pred2;
    
    
    subplot(2,num_outs,i);hold all;
    plot(out_bl); plot(pred_bl_p); plot(pred_bl_n);
    if do_corr
        c = corrcoef(out_bl,pred_bl_p); c1 = c(1,2);
        c = corrcoef(out_bl,pred_bl_n); c2 = c(1,2);
    else
        c1 = CalculateR2(out_bl,pred_bl_p);
        c2 = CalculateR2(out_bl,pred_bl_n);
    end
    
    all_c_bl(i,1) = c1;
    all_c_bl(i,2) = c2;
    
    title(['P: ' num2str(c1,2) ', N: ' num2str(c2,2)],'FontSize',14);
    if i == 1, ylabel('Baseline','FontSize',14); end
    axis('tight'); set(gca,'Box','off','TickDir','out','FontSize',14);
    
    
    
    subplot(2,num_outs,i+num_outs);hold all;
    plot(out_bl); plot(pred_ad_p); plot(pred_ad_n);
    if do_corr
        c = corrcoef(out_ad,pred_ad_p); c1 = c(1,2);
        c = corrcoef(out_ad,pred_ad_n); c2 = c(1,2);
    else
        c1 = CalculateR2(out_ad,pred_ad_p);
        c2 = CalculateR2(out_ad,pred_ad_n);
    end
    
    title(['P: ' num2str(c1,2) ', N: ' num2str(c2,2)],'FontSize',14);
    if i == 1, ylabel('Adaptation','FontSize',14); end
    axis('tight'); set(gca,'Box','off','TickDir','out','FontSize',14);
    
    all_c_ad(i,1) = c1;
    all_c_ad(i,2) = c2;
end

all_c_bl = reshape(all_c_bl_file,length(dates)*size(all_c_bl_file,1),2);
all_c_ad = reshape(all_c_ad_file,length(dates)*size(all_c_ad_file,1),2);


figure; hold all;
c=100*all_c_ad./all_c_bl;
m = mean(c,1);
s = [m-std(c,[],1); m+std(c,[],1)];
plot(m,'ko','LineWidth',3);
plot([1,2;1,2],s,'k-','LineWidth',2);
% plot(c','ko','LineWidth',3);
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[0 3]);

