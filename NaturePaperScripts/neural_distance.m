%% Prep datasets and get paths etc
clear; clc; close all;
dataSummary;

sessions      = { ...
    %     'Chewie','2016-09-09'; ... % VR
    %     'Chewie','2016-09-12'; ...
    %     'Chewie','2016-09-14'; ...
    %     'Chewie','2016-10-06'; ...
    %     'Mihili','2014-03-03'; ...
    %     'Mihili','2014-03-04'; ...
    %     'Mihili','2014-03-06'; ...
    'Chewie','2016-09-15'; ... % CF
    'Chewie','2016-10-05'; ...
    'Chewie','2016-10-07'; ...
    'Chewie','2016-10-11'; ...
%         'Mihili','2014-02-03'; ...
    'Mihili','2014-02-17'; ...
    'Mihili','2014-02-18'; ...
%     'Mihili','2014-03-07'; ...
%         'MrT','2013-08-19'; ... % CF
%         'MrT','2013-08-21'; ...
%         'MrT','2013-08-23'; ...
    %     'MrT','2013-09-03'; ... %VR
    %     'MrT','2013-09-05'; ...
    %     'MrT','2013-09-09'; ...
    };

% Session parameters
monkeys      = unique(sessions(:,1));
dates        = sessions(:,2);
tasks        = {'CO'};
perts        = {'FF'};


session_idx = getFileDBidx(filedb, ...
    {'Task',tasks,'Perturbation',perts,'Monkey',monkeys,'Date',dates}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});

idx_start = {'idx_target_on', 0};
idx_end   = {'idx_trial_end', -10};

badtrial_params = struct();%...
%     'ranges', {{'idx_go_cue','idx_movement_on',[5 40]; ...
%     'idx_peak_speed','idx_trial_end',[60 100]}});

badneuron_params_m1 = struct( ...
    'arrays','M1', ...
    'min_fr',1, ...
    'do_shunt_check',1, ...
    'use_trials',{{'epoch','BL'}});

badneuron_params_pmd = struct( ...
    'arrays','PMd', ...
    'min_fr',1, ...
    'do_shunt_check',1, ...
    'use_trials',{{'epoch',{'BL'}}});

func_calls = { ...
    {@getTDidx,'result','R'}, ...
    {@removeBadTrials,badtrial_params}, ...
    {@removeBadNeurons,badneuron_params_m1}, ...
    {@removeBadNeurons,badneuron_params_pmd}, ...
    @sqrtTransform, ...
    {@smoothSignals,struct('signals',{{'M1_spikes','PMd_spikes'}},'calc_fr',true,'kernel_SD',0.1)}, ...
    {@trimTD,idx_start,idx_end}, ...
    {@binTD,5}};

% load it
fname = cell(1,length(session_idx));
for iFile = 1:length(session_idx)
    file = session_idx(iFile);
    fname{iFile} = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
end
[trial_data,params] = loadTDfiles(fname, func_calls{:});


%%
close all;

which_signals = {'M1_pca','PMd_pca','PMdM1_potent','PMdM1_null'};
dims = 3;
bins_before = 1;
bins_after = 1;

figure;
for sig = 1:length(which_signals)
    which_signal= which_signals{sig};
    
    file_dist = cell(1,length(session_idx));
    for iFile = 1:length(session_idx)
        file = session_idx(iFile);
        
        [~,td] = getTDidx(trial_data,'date',datestr(filedb.Date{file},'mm-dd-yyyy'));
        % Find the potent and null spaces
        pn_params = struct( ...
            'in_signals','PMd_spikes', ...
            'out_signals','M1_spikes', ...
            'in_dims',16, ...
            'out_dims',8, ...
            'use_trials',{{'epoch',{'BL'}}}, ...
            'sqrt_transform',false, ...
            'do_smoothing',false);
        [td,pn_info] = getPotentSpace(td,pn_params);
        
        
        [~,temp_td] = getTDidx(td,'epoch','BL','range',[1 100]);
        u = unique([temp_td.target_direction]);
        bl_state = cell(1,length(u));
        for i = 1:length(u)
            idx = getTDidx(temp_td,'target_direction',u(i));
            temp_bl = zeros(length(idx),dims);
            for j = 1:length(idx)
                trial = idx(j);
                temp = temp_td(trial).(which_signal);
                temp_bl(j,:) = mean(temp(temp_td(trial).idx_go_cue-bins_before:temp_td(trial).idx_go_cue+bins_after,1:dims),1);
            end
            bl_state{i} = temp_bl;
        end
        
        idx = getTDidx(td,'epoch','BL');
        ad_state = zeros(length(idx),dims);
        ad_targs = zeros(length(idx),1);
        d = zeros(length(idx),1);
        for i = 1:length(idx)
            trial = idx(i);
            temp = td(trial).(which_signal);
            ad_state(i,:) = mean(temp(td(trial).idx_go_cue-bins_before:td(trial).idx_go_cue+bins_after,1:dims),1);
            ad_targs(i) = td(trial).target_direction;
            
            targ_id = find(u==ad_targs(i));
            
            temp = bl_state{targ_id};
            d(i) = mahal(ad_state(i,:),temp(:,1:dims));
        end
        
%         d = (mean(d(ceil(length(d)/2):end))-d)./mean(d(ceil(length(d)/2):end));
        
        file_dist{iFile} = d;
    end
    
    subplot(2,2,sig); hold all;
    for i = 1:length(file_dist)
        plot(moving_average(file_dist{i},30))
    end
    set(gca,'Box','off','TickDir','out','FontSize',14);
    title(which_signal);
    
    if sig == 1 || sig == 3
        ylabel('% Change in Mahal Dist.');
    end
    if sig == 3 || sig == 4
        xlabel('Curl Field Trials');
    end
end
