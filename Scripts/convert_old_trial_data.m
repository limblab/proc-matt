load('F:\TrialDataFiles\Chewie_CO_FF_2013-12-04.mat')

trial_data = trial_data(cellfun(@(x) ~isempty(x), {trial_data.idx_peak_speed}));

for i = 1:length(trial_data)
    trial_data(i).idx_trial_end = trial_data(i).idx_reward;
    trial_data(i).M1_spikes = full(trial_data(i).M1_spikes)';
    trial_data(i).PMd_spikes = full(trial_data(i).PMd_spikes)';
end



bad_idx = cellfun(@(x) any(sqrt(x(:,1).^2+x(:,2).^2) > 50),{trial_data.vel});
trial_data(bad_idx) = [];
disp(['Pruning ' num2str(sum(bad_idx)) ' trials with crazy velocities...']);

disp('Checking sorted units...');
trial_data = getCommonUnits(trial_data);
save('F:\TrialDataFiles\Chewie_CO_FF_2013-12-04.mat','trial_data');
