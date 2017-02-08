%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the root directory and paths
% defines root directories for all my computers
[~,hostname]= system('hostname');

switch strtrim(hostname)
    case 'FSMC17SW0D9GTF1' % my 2016 Macbook Pro
        rootDir = '/Users/mattperich/Data/';
    case 'FSM6YVJWR1' % my lab desktop
        rootDir = 'F:\';
    otherwise
        error('Computer not recognized. Could not pull root directory.')
end

cerebusDataDir = 'CerebusData';
CDSDir = 'CDS';
TDDir = 'TrialDataFiles';
resultsDir = 'results';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%