function [istuned, master_sg] = excludeCells(params,data,tuning,tracking,useArray)
% This function will check all of the cells against various exclusion
% criteria. For each cell:
%
% Cell-related criteria:
%   1) Is there a high SNR for waveforms
%           Look at waveforms for each cell
%   2) Is there a low percentage of short ISIs
%           Look at ISI of all spikes
%   3) Is there a minimum task-related average firing rate?
%           Look at FR used for tuning and take average
%   4) Is it the same cell in each epoch?
%
% Tuning-related criteria:
%   1) Is CI on each cell less than some level?
%   2) Is there an agreeable R2 for cosine fits?
%
%
% Data input is from baseline only. I assume that if it meets criteria
% there then it will for the rest, or else it won't pass the "same neuron"
% test.

tracking = tracking.(useArray);
data = data.(useArray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load some of the analysis parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ciSig = params.units.ciSignificance;
confLevel = params.classes.classConfidenceLevel;
numIters = params.tuning.numberBootIterations;
isiPercent = params.units.isiThreshold;
waveSNR = params.units.waveformSNR;
classifierBlocks = params.classes.classifierBlocks;
% a couple of parameters depend on the array
r2Min = params.units.([lower(useArray) '_r2Min']);
minFR = params.units.([lower(useArray) '_minFR']);
if isfield(params.units,[lower(useArray) '_isiPercent'])
    isiThresh = params.units.([lower(useArray) '_isiPercent']);
else
    isiThresh = params.units.isiPercent;
end

if isfield(tuning(1),'sg')
    all_sg = {tuning.sg};
else
    all_sg = {tuning.unit_guide};
end

tracking_chan = tracking{1}.chan;

badUnits = checkUnitGuides(all_sg);

% remove that index from each
[master_sg, ~] = setdiff(all_sg{1},badUnits,'rows');

for iBlock = 1:length(all_sg)
    [~, all_idx{iBlock}] = setdiff(all_sg{iBlock},badUnits,'rows');
end

% one element for each criteria. must meet criteria in all three task epochs
%   1) Waveform SNR
%   2) ISI Percentage
%   3) FR threshold
%   4) Neuron Tracking
%   5) PD CI
%   6) Cosine R2
istuned = zeros(size(master_sg,1),6);

units = data.units(all_idx{1});

%% Check that SNR of waveforms is above some threshold
% Compare peak to std of first bin, where presumably no cell is active
for unit = 1:size(master_sg,1)
    wf = units(unit).wf;
    %istuned(unit,1) = max(rms(wf'))/mean(std(double(wf(1:5,:)'))) >= waveSNR;
    sig = max(mean(wf'))-min(mean(wf'));
    istuned(unit,1) = sig / (2*mean(std(double(wf')))) >= waveSNR;
end


%% Check that ISI of neuron is above some threshold
for unit = 1:size(master_sg,1)
    isi = diff(units(unit).ts);
    istuned(unit,2) = sum(isi < isiThresh/1000)/length(isi) < isiPercent/100;
end


%% Check that neuron meets firing rate criterion in each epoch
if isfield(tuning(1),'fr')
    temp = all_idx{1};
    for unit = 1:size(master_sg,1)
        istuned(unit,3) = mean(tuning(1).fr(:,temp(unit)),1) >= minFR;
    end
else
    istuned(:,3) = ones(size(istuned(:,3)));
end

%% Check that same neuron is in each epoch
for unit = 1:size(master_sg,1)
    % Look at tracking output and determine if any are different
    idx = tracking_chan(:,1)==master_sg(unit,1)+.1*master_sg(unit,2);
    istuned(unit,4) = ~any(diff(tracking_chan(idx,:)));
end

%% don't do this if the classifierBlocks input is empty
if ~isempty(classifierBlocks)
    if isfield(tuning(1),'pds') % this means it is regression/glm/etc
        %% Check confidence in PD estimates
        all_pds = {tuning.pds};
        for unit = 1:size(master_sg,1)
            sig = zeros(size(all_pds));
            for iBlock = 1:length(all_pds)
                temp = all_pds{iBlock};
                temp = temp(all_idx{iBlock},:);
                sig(iBlock) = checkTuningCISignificance(temp(unit,:),ciSig,true);
            end
            istuned(unit,5) = all(sig(classifierBlocks));
        end
    elseif isfield(tuning(1),'mfr') % this means it is nonparametric
        % basically, here I want to check if any one bin is significantly
        % different from any other. I will call this a coarse estimate of
        % "tuning"
        nTargs = length(tuning(1).utheta);
        for unit = 1:size(master_sg,1)
            blockDiffs = zeros(1,length(tuning));
            for iBlock = 1:length(tuning)
                cil = tuning(iBlock).cil;
                cih = tuning(iBlock).cih;
                
                targDiffs = [];
                for iTarg = 1:nTargs-1
                    % get confidence intervals for
                    for iTarg2 = iTarg+1:nTargs
                        % do the confidence intervals overlap?
                        ci1 = [cil(unit,iTarg) cih(unit,iTarg)];
                        ci2 = [cil(unit,iTarg2) cih(unit,iTarg2)];
                        overlap = range_intersection(ci1,ci2);
                        
                        % build matrix showing how they differ
                        if isempty(overlap)
                            targDiffs = [targDiffs; 1];
                        else
                            targDiffs = [targDiffs; 0];
                        end
                    end
                end
                blockDiffs(iBlock) = sum(targDiffs) > 2;
            end
            %istuned(unit,5) = all(blockDiffs);
            p = anova1(tuning(iBlock).boot_fr{unit},[],'off');
            istuned(unit,5) = p < 0.05;
        end
    end
    %% Check that r-squared of fit is okay
    if isfield(tuning(1),'r_squared') && ~isempty(tuning(1).r_squared)
        all_rs = {tuning.r_squared};
        all_rs_ci = cell(size(all_rs));
        
        for iBlock = 1:length(all_rs)
            temp = all_rs{iBlock};
            temp = sort(temp(all_idx{iBlock},:),2);
            
            all_rs{iBlock} = temp;
            
            % get 95% CI for each
            all_rs_ci{iBlock} = [temp(:,ceil(numIters - (confLevel/2)*numIters)), temp(:,floor((confLevel/2)*numIters))];
        end
        
        for unit = 1:size(master_sg,1)
            sig = zeros(size(all_rs));
            for iBlock = 1:length(all_rs)
                % also only consider cells that are described by cosines
                %   have bootstrapped r2... see if 95% CI is > threshold?
                temp = all_rs_ci{iBlock};
                sig(iBlock) = temp(unit,1) > r2Min;
            end
            
            % check significance
            % only consider cells that are tuned in all epochs
            %   first column is CI bound, second is r-squared
            istuned(unit,6) = all(sig(classifierBlocks));
        end
    else
        istuned(:,6) = ones(size(istuned(:,6)));
    end
    
end

