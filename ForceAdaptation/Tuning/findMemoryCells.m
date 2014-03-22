function classes = findMemoryCells(expParamFile, outDir, paramSetName, compMethod)
% FINDMEMORYCELLS  Compares tuning of cells to classify their behavior
%
%   This function uses the tuning made by fitTuningCurves to classify them
% based on their tuning behavior before, during, and after learning. There
% are the following types:
%       Kinematic (AAA, id=1) : no change with perturbation
%       Dynamic   (ABA, id=2)   : changes in response to the perturbation
%       Memory I  (ABB, id=3)  : changes with perturbation and hold tuning
%       Memory II (AAB, id=4) : changes when perturbation is removed
%       Other     (ABC, id=5)     : Different tuning in every epoch
%
% INPUTS:
%   expParamFile: (string) path to file containing experimental parameters
%   outDir: (string) directory for output
%   paramSetName: (string) Name of the parameter file
%   compMethod: (string) (optional) 'diff' or 'overlap', method of sig
%
% OUTPUTS:
%   classes: (struct) output with field for each array and tuning period
%
% NOTES:
%   - Assumes the Baseline->Adaptation->Washout files exist and uses all 3
%   - This function requires several bits of pre-processing
%       1) Create a data struct from the Cerebus files (makeDataStruct)
%       2) Fit tuning for neurons, regression and nonparametric recommended (fitTuningCurves)
%   - This function will automatically write the struct to a file, too
%   - See "experimental_parameters_doc.m" for documentation on expParamFile
%   - Analysis parameters file must exist (see "analysis_parameters_doc.m")

if nargin < 4
    compMethod = 'diff';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load some of the experimental parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = parseExpParams(expParamFile);
useDate = params.date{1};
taskType = params.task{1};
adaptType = params.adaptation_type{1};
clear params

dataPath = fullfile(outDir,useDate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load some of the analysis parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paramFile = fullfile(dataPath, paramSetName, [ useDate '_' paramSetName '_tuning_parameters.dat']);
params = parseExpParams(paramFile);
tuningPeriods = params.tuning_periods;
tuningMethods = params.tuning_methods;
clear params;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveFile = fullfile(dataPath,paramSetName,[taskType '_' adaptType '_classes_' useDate '.mat']);

disp('Loading data to classify cells...')
tuningFile = fullfile(dataPath,paramSetName,[taskType '_' adaptType '_tuning_' useDate '.mat']);
load(tuningFile,'BL');
load(tuningFile,'AD');
load(tuningFile,'WO');

arrays = BL.meta.arrays;

for iArray = 1:length(arrays)
    useArray = arrays{iArray};
    for iMethod = 1:length(tuningMethods)
        % nonparametric tuning requires a different method for comparison
        % for regression or GLM, loop along the periods
        for iPeriod = 1:length(tuningPeriods)
            % only glm can use the full file tuning
            if strcmpi(tuningPeriods{iPeriod},'file') && ~strcmpi(tuningMethods{iMethod},'glm')
                warning(['File tuning not supported for ' tuningMethods{iMethod} ' method...']);
            elseif strcmpi(tuningMethods{iMethod},'nonparametric')
                warning(['Classification not supported for ' tuningMethods{iMethod} ' method...']);
            else
                for iBlock = 1:length(AD.(useArray).(tuningMethods{iMethod}).(tuningPeriods{iPeriod}))
                    blt = BL.(useArray).(tuningMethods{iMethod}).(tuningPeriods{iPeriod});
                    adt = AD.(useArray).(tuningMethods{iMethod}).(tuningPeriods{iPeriod});
                    wot = WO.(useArray).(tuningMethods{iMethod}).(tuningPeriods{iPeriod});
                    [cellClass,sg] = classifyCells(blt,adt(iBlock),wot,BL.meta,tuningMethods{iMethod},compMethod,paramSetName);
                    
                    % get cells that are significantly tuned in all epochs
                    %   first column is PDs, second is MDs
                    tunedCells = find(cellClass(:,1)~=-1);
                    disp(['There are ' num2str(length(tunedCells)) ' cells tuned in all epochs...']);
                    
                    s(iBlock).classes = cellClass;
                    s(iBlock).sg = sg;
                    s(iBlock).tuned_cells = tunedCells;
                end
                classes.(useArray).(tuningMethods{iMethod}).(tuningPeriods{iPeriod}) = s;
            end
        end
    end
end

% save the new file with classification info
disp(['Saving data to ' saveFile]);
save(saveFile,'-struct','classes');
