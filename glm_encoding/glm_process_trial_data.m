function [trial_data, params] = glm_process_trial_data(trial_data,params)

fn = fieldnames(params);
for i = 1:length(fn)
    eval([fn{i} ' = params.' fn{i} ';']);
end, clear i good_cells cov_array pred_array pert epochs pca_w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate speed and target direction if desired
if any(ismember(kin_signals,'speed'))
    for trial = 1:length(trial_data)
        trial_data(trial).speed = hypot(trial_data(trial).vel(:,1), ...
            trial_data(trial).vel(:,2));
    end
end
if any(ismember(kin_signals,'targ'))
    trial_data = getTargetDirection(trial_data);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter out neurons
trial_data = getCommonUnits(trial_data);

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
        idx = trial_data(iTrial).(train_start_idx{1})+train_start_idx{2}:trial_data(iTrial).(train_end_idx{1})+train_end_idx{2};
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
if do_rcb
    % Convolve with basis functions before
    disp('Convolving with raised cosines...');
    rcb_which_vars = cell(1,length(arrays));
    for i = 1:length(arrays)
        rcb_which_vars{i} = [arrays{i} '_spikes'];
    end
    rcb_which_vars = [rcb_which_vars, kin_signals];
    params.rcb_which_vars = rcb_which_vars;
    params.rcb_n = params.unit_lags;
    tic;
    trial_data = convBasisFunc(trial_data,rcb_which_vars,params);
    toc;
    
    for i = 1:length(kin_signals)
        % shift kinematics by kin_lags backwards
        for j = 1:length(trial_data)
            temp = trial_data(j).(kin_signals{i});
            trial_data(j).(kin_signals{i}) = [temp(kin_lags+1:end,:); NaN(kin_lags,size(temp,2))];
            temp = trial_data(j).([kin_signals{i} '_shift']);
            trial_data(j).([kin_signals{i} '_shift']) = [temp(kin_lags+1:end,:); NaN(kin_lags,size(temp,2))];
        end
    end
elseif ~isempty(kin_signals) || do_all_history || do_self_history
    error('FIX KINEMATIC SHIFT');
    if unit_lags > 0
        % Duplicate and shift
        build_inputs = cell(1,2*(length(kin_signals) + length(arrays)));
        for i = 1:length(arrays)
            build_inputs{(i-1)*2+1} = [arrays{i} '_spikes'];
            build_inputs{(i-1)*2+2} = unit_lags;
        end
        if ~isempty(kin_signals)
            for i = 1:length(kin_signals)
                build_inputs{2*length(arrays) + (i-1)*2+1} = kin_signals{i};
                build_inputs{2*length(arrays) + (i-1)*2+2} = max(kin_lags);
            end
        end
        trial_data = dupeAndShift(trial_data,build_inputs);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
