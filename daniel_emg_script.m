load('/Users/mattperich/Data/Han_20170206_CObumpcurl_baseline_001_area2EMG_CDS.mat');
td1 = parseFileByTrial(cds,struct('array_alias',{{'LeftS1Area2','S1'}},'exclude_units',255,'event_list',{{'goCueTime';'bumpTime';'bumpDir'}}));
load('/Users/mattperich/Data/Han_20170206_CObumpcurl_baseline_002_area2EMG_CDS.mat');
[td2,td_params] = parseFileByTrial(cds,struct('array_alias',{{'LeftS1Area2','S1'}},'exclude_units',255,'event_list',{{'goCueTime';'bumpTime';'bumpDir'}}));
trial_data = [td1 td2];
save('/Users/mattperich/Data/han_td.mat','trial_data','td_params');

%%
clear;
clc;
close all;

% Load data
load('/Users/mattperich/Data/han_td.mat','trial_data');

%%
utheta = unique([trial_data.target_direction]);
td = truncateAndBin(trial_data,{'idx_goCueTime',-20},{'idx_goCueTime',20});
% td = trialAverage(td,'target_direction');
[~,td] = getTDidx(td,'target_direction',utheta(1));

idx = 4;
figure;
hold all;
for i = 1:length(td)
    plot(td(i).emg(:,idx));
end

