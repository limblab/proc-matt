close all;
clear;
clc;

dataSummary;

filenames = { ...
    'Chewie_CO_FF_2013-10-22', ...
    'Chewie_CO_FF_2013-10-23', ...
    'Chewie_CO_FF_2013-10-31', ...
    'Chewie_CO_FF_2013-11-01', ...
    'Chewie_CO_FF_2013-12-03', ...
%     'Chewie_CO_FF_2013-12-04', ...
    'Chewie_CO_FF_2015-07-01', ...
    'Chewie_CO_FF_2015-07-03', ...
        'Chewie_CO_FF_2015-06-29', ...
    'Chewie_CO_FF_2015-06-30', ...
    'Chewie_CO_FF_2015-07-06', ...
%     'Chewie_CO_FF_2015-07-07', ...
%     'Chewie_CO_FF_2015-07-08', ...
%     'Mihili_CO_FF_2014-02-03', ...
%     'Mihili_CO_FF_2014-02-17', ...
%     'Mihili_CO_FF_2014-02-18', ...
%     'Mihili_CO_FF_2014-03-07', ...
%     'Mihili_CO_FF_2015-06-11', ...
%     'Mihili_CO_FF_2015-06-17', ...
%     'Mihili_CO_FF_2015-06-10', ...
%     'Mihili_CO_FF_2015-06-15', ...
%     'Mihili_CO_FF_2015-06-16', ...
    };

dataDir = fullfile(rootDir,TDDir);

total_counts = zeros(length(filenames),1);

trial_counts = zeros(length(filenames),3);
for iFile = 1:length(filenames)
    load(fullfile(dataDir,[filenames{iFile} '.mat']),'trial_data');
    total_counts(iFile) = length(trial_data);
    if isfield(trial_data,'result'), trial_data = trial_data(getTDidx(trial_data,'result',{'A','I','F'})); end
    
    trial_counts(iFile,1) = sum(getTDidx(trial_data,'epoch','BL'));
    trial_counts(iFile,2) = sum(getTDidx(trial_data,'epoch','AD'));
    trial_counts(iFile,3) = sum(getTDidx(trial_data,'epoch','WO'));
end

max(sum(trial_counts,2)./total_counts)