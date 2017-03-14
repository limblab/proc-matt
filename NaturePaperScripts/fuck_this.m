function trial_data = fuck_this(trial_data)

load('/Users/mattperich/Data/results/trainad_null_bl/FF-PMd-M1_Chewie_CO_FF_2016-10-07.mat','params');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% populate trial_data struct with PCA scores and null/potent projections
    for trial = 1:length(trial_data)
        % get signals to recreate input PCA
        data = get_vars(trial_data(trial),{'PMd_spikes',params.good_cells{2}});
        trial_data(trial).PMd_pca = sqrt(data) * params.pca_w;
        
            % now do null/potent
            data = trial_data(trial).PMd_pca;
            data = data(:,params.pca_dims.PMd);
            trial_data(trial).PMdM1_potent = data * params.V_potent;
            trial_data(trial).PMdM1_null = data * params.V_null;
    end
