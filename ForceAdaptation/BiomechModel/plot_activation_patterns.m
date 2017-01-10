% make a Lillicrap 2013 esque plot of muscle activations as function of
% time and direction
clear
clc;
close all;

load('F:\TrialDataFiles\Mihili_CO_FF_2014-02-17.mat')
load('F:\TrialDataFiles\biomech_model\Mihili_CO_FF_2014-02-17_sim.mat')

if isfield(trial_data(1),'result'),
    trial_data = trial_data(strcmpi({trial_data.result},'R'));
end


%     idx = find(getTDidx(trial_data,'epoch','AD'));
%     trial_data = trial_data(idx);

% good_idx = cellfun(@(x) sum(x~=1),{sim_data.muscle_flags}) == 0;
% sim_data = sim_data(good_idx);
% trial_data = trial_data(good_idx);

%%
utheta = unique([trial_data.target_direction]);
muscle_names = {'Shoulder-Flex','Shoulder-Ext','Elbow-Flex','Elbow-Ext','BiArt-Flex','BiArt-Ext'};

plot_order = [1,4,2,5,3,6];

figure;
for m = 1:6
    [targs,targs2] = deal(cell(1,length(utheta)));
    for i = 1:length(utheta)
        [targs{i},targs2{i}] = deal([]);
    end
    
    s = sim_data(getTDidx(trial_data,'epoch','BL'));
    t = trial_data(getTDidx(trial_data,'epoch','BL'));
    
    for trial = 1:length(s)
        % get indices of go cue and trial end
        idx = t(trial).idx_movement_on:t(trial).idx_trial_end-50;
        
        t_idx = find(utheta == t(trial).target_direction);
        
        
        temp = s(trial).muscles(:,m);
        % interpolate to 100 datapoints
        targs{t_idx} = [targs{t_idx}; interp1(1:length(idx),temp(idx),linspace(1,length(idx),100))];
    end
    
    idx = find(getTDidx(trial_data,'epoch','AD'));
    idx = idx(floor(0.5*length(idx)):end);
    s = sim_data(idx);
    t = trial_data(idx);
    
    for trial = 1:length(s)
        % get indices of go cue and trial end
        idx = t(trial).idx_movement_on:t(trial).idx_trial_end-50;
        
        t_idx = find(utheta == t(trial).target_direction);
        
        
        temp = s(trial).muscles(:,m);
        % interpolate to 100 datapoints
        targs2{t_idx} = [targs2{t_idx}; interp1(1:length(idx),temp(idx),linspace(1,length(idx),100))];
    end
    
    
    subplot(2,3,plot_order(m));
    
    ma1 = cell2mat(cellfun(@(x) mean(x,1),targs,'Uni',0)');
%     ma1 = ma1./repmat(max(ma1,[],2),1,size(ma1,2));
    ma2 = cell2mat(cellfun(@(x) mean(x,1),targs2,'Uni',0)');
%     ma2 = ma2./repmat(max(ma2,[],2),1,size(ma1,2));
    
    imagesc(1:100,utheta*180/pi,ma2);
    axis('tight');
    set(gca,'Box','off','TickDir','out','FontSize',14);
    if plot_order(m)>3
        xlabel('Normalized Time','FontSize',14);
    else
        set(gca,'XTickLabel',[]);
    end
    if plot_order(m)==1 || plot_order(m)==4
        ylabel('Target Direction','FontSize',14);
    else
        set(gca,'YTickLabel',[]);
    end
    title(muscle_names{m},'FontSize',14);
end



%% this one is neural

% Theories of optimal control suggest you don't need any explicit
% representation. During learning, you can compute a new optimal motor plan
% without changing how motor cortex drives the limb. We present results
% showing that during learning, M1 maintains a fixed relationship with
% dynamics. We present a basic model that can replicate the tuning
% properties of our recorded neurons while maintaining fixed relationships
% with dynamical variables.

utheta = unique([trial_data.target_direction]);
muscle_names = {'Shoulder-Flex','Shoulder-Ext','Elbow-Flex','Elbow-Ext','BiArt-Flex','BiArt-Ext'};

bin_size            = 0.01;
kernel_SD           = 3*bin_size;

% kernel half length is 3·SD out
kernel_hl               = ceil( 3 * kernel_SD / (bin_size) );
% create the kernel --it will have length 2*kernel_hl+1
kernel                  = normpdf( -kernel_hl*(bin_size) : ...
    bin_size : kernel_hl*(bin_size), ...
    0, kernel_SD );
% compute normalization factor --this factor depends on the number of taps
% actually used

if 0
    n_neurons = size(trial_data(1).M1_spikes,1);
else
    n_neurons = size(sim_data(1).muscle_neurons,1);
end

for unit = 1:n_neurons
    [targs,targs2] = deal(cell(1,length(utheta)));
    for i = 1:length(utheta)
        [targs{i},targs2{i}] = deal([]);
    end
    
    s = sim_data(getTDidx(trial_data,'epoch','BL'));
    t = trial_data(getTDidx(trial_data,'epoch','BL'));
    
    for trial = 1:length(s)
        % get indices of go cue and trial end
        idx = t(trial).idx_movement_on-20:t(trial).idx_trial_end-50;
        if 0
            nm = conv(kernel,ones(1,size(t(trial).M1_spikes,1)))';
        else
            nm = conv(kernel,ones(1,size(s(trial).muscle_neurons,1)))';
        end
        
        t_idx = find(utheta == t(trial).target_direction);
        
        if 0
            temp = t(trial).M1_spikes(:,unit);
        else
            temp = s(trial).muscle_neurons(:,unit);
        end
        temp = conv(kernel,temp) ./ nm;
        temp = temp(kernel_hl+1:end-kernel_hl);
        
        
        % interpolate to 100 datapoints
        targs{t_idx} = [targs{t_idx}; interp1(1:length(idx),temp(idx),linspace(1,length(idx),1000))];
    end
    
    figure;
    imagesc(1:1000,utheta*180/pi,cell2mat(cellfun(@(x) mean(x,1),targs,'Uni',0)'));
    pause;
    close all;
end

