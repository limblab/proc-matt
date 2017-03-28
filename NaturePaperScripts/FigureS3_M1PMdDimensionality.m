%% Prep datasets and get paths etc
clear; clc; close all;
dataSummary;
trial_params;

sessions = { ...
    'Chewie','2016-09-09'; ... % VR
    'Chewie','2016-09-12'; ...
    'Chewie','2016-09-14'; ...
    'Chewie','2016-10-06'; ...
    'Mihili','2014-03-03'; ...
    'Mihili','2014-03-04'; ...
    'Mihili','2014-03-06'; ...
    'Chewie','2016-09-15'; ... % CF
    'Chewie','2016-10-05'; ...
    'Chewie','2016-10-07'; ...
    'Chewie','2016-10-11'; ...
    'Mihili','2014-02-03'; ...
    'Mihili','2014-02-17'; ...
    'Mihili','2014-02-18'; ...
    'Mihili','2014-03-07'};


% Session parameters
monkeys = unique(sessions(:,1));
dates = unique(sessions(:,2));
tasks = {'CO'};
perts = {'FF','VR'};

session_idx = getFileDBidx(filedb, ...
    {'Date',dates,'Task',tasks,'Perturbation',perts,'Monkey',monkeys});%, ...
    %{'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    %'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});

% Load data
fnames = cell(1,length(session_idx));
for iFile = 1:length(session_idx)
    file = session_idx(iFile);
    fnames{iFile} = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
end

func_calls = [trial_func_calls { ...
    {@removeBadNeurons,badneuron_params}, ...
    {@removeBadTrials,badtrial_params}, ...
    {@getTDidx,'epoch','BL'}, ...
    @sqrtTransform, ...
    {@smoothSignals,struct('signals',{{'M1_spikes','PMd_spikes'}},'kernel_SD',0.1)}, ...
    }];

[trial_data,params] = loadTDfiles(fnames, func_calls{:});


%%
figure;
arrays = {'M1','PMd'};
max_dim = 15;

all_dims = cell(1,length(arrays));
for a = 1:length(arrays)
    
    array = arrays{a};
    
    dims = zeros(1,length(session_idx));
    for i = 1:length(session_idx)
        disp(['Working on ' array ', file ' num2str(i) ' of ' num2str(length(session_idx))]);
        file = session_idx(i);
        
        [~,td] = getTDidx(trial_data, ...
            'date',datestr(filedb.Date{file},'mm-dd-yyyy'), ...
            'monkey',filedb.Monkey{file}, ...
            'task',filedb.Task{file},...
            'perturbation',filedb.Perturbation{file});
        
        if strcmpi(array,'M1')
            td = trimTD(td,{'idx_movement_on',-20},{'idx_movement_on',50});
        else
            td = appendTDs( ...
            trimTD(td,{'idx_target_on',0},{'idx_target_on',50}), ...
            trimTD(td,{'idx_go_cue',0},{'idx_go_cue',20}), ...
            trimTD(td,{'idx_movement_on',-10},{'idx_movement_on',50}));
        end
        
        
        for j = 1:length(td)
            td(j).target_direction = bin_angles(td(j).target_direction,2*pi/8);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ESTIMATE DIMENSIONALITY
        dims(i) = estimateDimensionality(td,struct('signals',[array '_spikes'],'condition','target_direction'));
        clear td_temp
    end
    all_dims{a} = dims;
    
    subplot(length(arrays),1,a);
    hist(dims,1:max_dim);
    axis('tight');
    set(gca,'Box','off','TickDir','out','FontSize',14);
    ylabel('Counts','FontSize',14);
    title([array ' (N = ' num2str(length(dims)) ')'],'FontSize',16);
end

xlabel('Dimensionality','FontSize',14);