%% Prep datasets and get paths etc
clear; clc; close all;
dataSummary;
trial_params;

% Session parameters
monkeys = {'Chewie','Mihili'};
tasks = {'CO'};
perts = {'FF'};

session_idx = getFileDBidx(filedb, ...
    {'Task',tasks,'Perturbation',perts,'Monkey',monkeys}, ...
    {'~(ismember(filedb.Monkey,''Mihili'') & datenum(filedb.Date) > datenum(''2015-01-01''))', ...
    'cellfun(@(x) all(ismember({''M1'',''PMd''},x)),filedb.Arrays)'});

% Load data
fnames = cell(1,length(session_idx));
for iFile = 1:length(session_idx)
    file = session_idx(iFile);
    fnames{iFile} = fullfile(rootDir,TDDir,[filedb.Monkey{file} '_' filedb.Task{file} '_' filedb.Perturbation{file} '_' filedb.Date{file} '.mat']);
end

func_calls = [trial_func_calls { ...
    {@removeBadNeurons,badneuron_params}, ...
    {@getTDidx,'epoch','BL'}, ...
    @sqrtTransform, ...
    {@smoothSignals,struct('signals',{{'M1_spikes','PMd_spikes'}},'kernel_SD',pn_kernel_SD)}, ...
    {@trimTD,{'idx_go_cue',-50},{'idx_go_cue',50}}}];

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