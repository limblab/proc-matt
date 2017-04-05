

dataSummary;
dataDir = fullfile(rootDir,TDDir);

filenames = { ...
    'Chewie_CO_FF_2013-10-22', ...
    'Chewie_CO_FF_2013-10-23', ...
    'Chewie_CO_FF_2013-10-31', ...
    'Chewie_CO_FF_2013-11-01', ...
    'Chewie_CO_FF_2013-12-03', ...
    'Chewie_CO_FF_2013-12-04', ...
    'Chewie_CO_FF_2015-07-01', ...
    'Chewie_CO_FF_2015-07-03', ...
    'Chewie_CO_FF_2015-06-29', ...
    'Chewie_CO_FF_2015-06-30', ...
    'Chewie_CO_FF_2015-07-06', ...
    'Chewie_CO_FF_2015-07-07', ...
    'Chewie_CO_FF_2015-07-08', ...
    'Mihili_CO_FF_2014-02-03', ...
    'Mihili_CO_FF_2014-02-17', ...
    'Mihili_CO_FF_2014-02-18', ...
    'Mihili_CO_FF_2014-03-07', ...
    'Mihili_CO_FF_2015-06-11', ...
    'Mihili_CO_FF_2015-06-17', ...
    'Mihili_CO_FF_2015-06-10', ...
    'Mihili_CO_FF_2015-06-15', ...
    'Mihili_CO_FF_2015-06-16', ...
    };

num_channels = 96;

chans_with_units = zeros(1,length(filenames));
chans_with_mult = zeros(1,length(filenames));

for iFile = 1:length(filenames)
    % load each trial_data file
    load(fullfile(dataDir,[filenames{iFile} '.mat']),'trial_data');
    
    % get unique electrodes
    e = unique(trial_data(1).M1_unit_guide(:,1));
    
    % get number of electrodes with more than one unit
    u = zeros(size(e));
    for i = 1:length(e)
        u(i) = sum(trial_data(1).M1_unit_guide(:,1) == e(i));
    end
    
    chans_with_units(iFile) = length(e);
    chans_with_mult(iFile) = sum(u > 1);
end

% chewie is first 13

