function tracking = trackNeuronsAcrossEpochs(params)
% TRACKNEURONS  Run empirical KS test to check for stability of neurons
%
%   This function will allow you to track cells across a session for the
% three epochs (baseline, adaptation, washout).
%
% INPUTS:
%   expParamFile: (string) path to file containing experimental parameters
%   outDir: (string) directory for output
%   criteria: (cell array of strings) for significance, 'isi' and/or 'wf'
%
% OUTPUTS:
%   tracking: (struct) results
%       Output may not be super intuitive to interpret so talk to Matt or
%         Brian for now until I can code up a description.
%
% NOTES:
%   -This function requires several bits of pre-processing
%       1) Create a data struct from the Cerebus files (makeDataStruct)
%   - This function will automatically write the struct to a file, too
%   - See "experimental_parameters_doc.m" for documentation on expParamFile
%   - Analysis parameters file must exist (see "analysis_parameters_doc.m")

procDirName = 'Processed';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load some of the experimental parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root_dir = params.outDir; % we want to load from the output directory of makeDataStruct
useDate = params.exp.date;
taskType = params.exp.task;
adaptType = params.exp.adaptation_type;
monkey = params.exp.monkey;
epochs = params.exp.epochs;

dataPath = fullfile(root_dir,procDirName,useDate);

if ~iscell(epochs)
    epochs = {epochs};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Loading data to track neurons...')
data = cell(1,length(epochs));
for iEpoch = 1:length(epochs)
    data{iEpoch} = loadResults(root_dir,{monkey, useDate, adaptType, taskType},'data',[],epochs{iEpoch});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveFile = fullfile(dataPath,[taskType '_' adaptType '_tracking_' useDate '.mat']);

tracking = trackNeurons(params,data);

% save the new file with classification info
disp(['Saving data to ' saveFile]);
save(saveFile,'-struct','tracking');

