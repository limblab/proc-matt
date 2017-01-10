clear; clc; close all;
dataSummary;

sessions = { ...
    'Chewie','2016-09-15'; ...
    'Chewie','2016-09-19'; ...
    'Chewie','2016-10-05'; ...
    'Chewie','2016-10-07'; ...
    'Mihili','2014-02-17'; ...
    'Mihili','2014-02-18'; ...
    'Mihili','2014-03-07'; ...
    'Mihili','2015-06-15'; ...
    'Mihili','2015-06-16'; ...
    'Mihili','2015-06-17'; ...
    'Mihili','2014-02-03'; ...
    };


sessions = { ...
    'Chewie','2016-09-09'; ...
    'Chewie','2016-09-12'; ...
    'Chewie','2016-09-14'; ...
    'Mihili','2014-03-03'; ...
    'Mihili','2014-03-04'; ...
    'Mihili','2014-03-06'; ...
    'Mihili','2015-06-23'; ...
    'Mihili','2015-06-25'; ...
    'Mihili','2015-06-26'; ...
    'MrT','2013-09-03'; ...
    'MrT','2013-09-05'; ...
    'MrT','2013-09-09'; ...
    };

sessions = sessions(3,:);

num_dims = 12;
rand_bl_runs = 5;
arrays = {'M1','PMd','M1PMd'};

monkeys = {'Chewie','Mihili'};
tasks = {'CO','RT'};
pert = 'VR';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
font_size = 14;
line_width = 2;
plot_colors = [0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    %     0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; subplot1(1,length(arrays));

use_date_idx = ismember(filedb.Monkey,monkeys) & ismember(filedb.Perturbation,pert) & ismember(filedb.Task,tasks);
use_date_idx = use_date_idx & ismember(filedb.Date,sessions(:,2));
use_files = find(use_date_idx);

for iFile = 1:length(use_files)
    filename = [filedb.Monkey{use_files(iFile)} '_' filedb.Task{use_files(iFile)} '_' filedb.Perturbation{use_files(iFile)} '_' filedb.Date{use_files(iFile)}];
    load(fullfile(rootDir,TDDir,[filename '.mat']));
    
    for iArray = 1:length(arrays)
        if strcmpi(arrays{iArray},'m1pmd')
            which_array = {'M1','PMd'};
        else
            which_array = arrays{iArray};
        end
        
        subplot1(iArray); hold all;
        
        % make a few random BL comparisons
        bl_idx = find(getTDidx(trial_data,'epoch','bl'));
        for i = 1:rand_bl_runs
            idx = randperm(length(bl_idx));
            w1 = getPCA(trial_data,struct('array',which_array,'bin_size',0.03,'trial_idx', bl_idx(idx(1:floor(length(idx)/2))) ));
            w2 = getPCA(trial_data,struct('array',which_array,'bin_size',0.03,'trial_idx', bl_idx(idx(floor(length(idx)/2)+1:end)) ));
            pa = principal_angles(w1(:,1:num_dims),w2(:,1:num_dims));
            plot(pa.*180/pi,'LineWidth',1,'Color',[0.75 0.75 0.75]);
        end
        
        w1 = getPCA(trial_data,struct('array',which_array,'bin_size',0.03,'trial_idx', bl_idx(1:floor(length(bl_idx)/2)) ));
        w2 = getPCA(trial_data,struct('array',which_array,'bin_size',0.03,'trial_idx', bl_idx(floor(length(bl_idx)/2)+1:end) ));
        pa = principal_angles(w1(:,1:num_dims),w2(:,1:num_dims));
        plot(pa.*180/pi,'LineWidth',line_width,'Color','k');
        
        [w_bl] = getPCA(trial_data,struct('array',which_array,'bin_size',0.03,'trial_idx',bl_idx));
        
        trial_idx = find(getTDidx(trial_data,'epoch','ad'));
        [w] = getPCA(trial_data,struct('array',which_array,'bin_size',0.03,'trial_idx',trial_idx(1:floor(length(bl_idx)/2))));
        pa = principal_angles(w2(:,1:num_dims),w(:,1:num_dims));
        plot(pa.*180/pi,'LineWidth',line_width,'Color',plot_colors(2,:));
        
        [w] = getPCA(trial_data,struct('array',which_array,'bin_size',0.03,'trial_idx',trial_idx(end-floor(length(bl_idx)/2):end)));
        pa = principal_angles(w2(:,1:num_dims),w(:,1:num_dims));
        plot(pa.*180/pi,'LineWidth',line_width,'Color',plot_colors(1,:));
        
        
        trial_idx = find(getTDidx(trial_data,'epoch','wo'));
        [w] = getPCA(trial_data,struct('array',which_array,'bin_size',0.03,'trial_idx',trial_idx(1:floor(length(bl_idx)/2))));
        pa = principal_angles(w2(:,1:num_dims),w(:,1:num_dims));
        plot(pa.*180/pi,'LineWidth',line_width,'Color',plot_colors(3,:));
        
        trial_idx = find(getTDidx(trial_data,'epoch','wo'));
        [w] = getPCA(trial_data,struct('array',which_array,'bin_size',0.03,'trial_idx',trial_idx(end-floor(length(bl_idx)/2):end)));
        pa = principal_angles(w2(:,1:num_dims),w(:,1:num_dims));
        plot(pa.*180/pi,'LineWidth',line_width,'Color',plot_colors(4,:));
        
        set(gca,'Box','off','TickDir','out','FontSize',font_size,'XLim',[1 num_dims],'YLim',[0 90]);
        if iArray == 1
            ylabel('Principal Angle','FontSize',font_size);
        end
        xlabel('Dimension','FontSize',font_size);
        title(arrays{iArray},'FontSize',font_size+2);
    end
end

% legend({[pert '-early'],[pert '-late'],'WO-early','WO-late'},'FontSize',14);

