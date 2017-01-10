function [trial_data,bad_trials] = filterTrials(trial_data,params)

% all of these are in number of bins [min,max]
reaction_time = [10 60];
angle_error = pi/2;
result_codes = {'R'};

if isfield(trial_data(1),'result')
    trial_data = trial_data(ismember({trial_data.result},result_codes));
end

[rt,ttt] = deal(zeros(1,length(trial_data)));
for trial = 1:length(trial_data)
    rt(trial) = trial_data(trial).idx_movement_on - trial_data(trial).idx_go_cue;
    ttt(trial) = trial_data(trial).idx_movement_on - trial_data(trial).idx_trial_end;
end

err = get_learning_metrics(trial_data,'angle');

bad_trials = rt < reaction_time(1) | rt > reaction_time(2) | abs(err') > angle_error;

trial_data = trial_data(~bad_trials);
