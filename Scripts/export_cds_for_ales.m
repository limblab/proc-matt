% load CDS
% load('F:\Chewie\CDS\2016-10-21\Chewie_CO_CS_BL_10212016_001.mat');

clear data meta trials lfp units;
targ_angs = [pi/2, pi/4, 0, -pi/4, -pi/2, -3*pi/4, pi, 3*pi/4];
%
meta = cds.meta;
meta = rmfield(meta,{'cdsVersion','processedTime','knownProblems','hasEmg','hasLfp','hasKinematics','hasForce','hasAnalog','hasUnits','hasTriggers','hasBumps','hasChaoticLoad','hasSorting','numSorted','percentStill','numDualUnits','numWellSorted','stillTime','aliasList','cdsName'});
data.meta = meta;

trials = cds.trials;
for trial = 1:size(cds.trials,1)
    if cds.trials.tgtID(trial)+1 <= 8
        trials.tgtDir(trial) = targ_angs(cds.trials.tgtID(trial)+1);
    else
        trials.tgtDir(trial) = NaN;
    end
end
trials.bumpTime = [];
trials.bumpID = [];
trials.bumpPhase = [];
trials.bumpDir = [];
data.trials = trials;

kin = struct( ...
    't',cds.kin.t, ...
    'x',cds.kin.x, ...
    'y',cds.kin.y, ...
    'vx',cds.kin.vx, ...
    'vy',cds.kin.vy, ...
    'ax',cds.kin.ax, ...
    'ay',cds.kin.ay);
data.kin = kin;

vn = cds.lfp.Properties.VariableNames;
for i = 1:length(vn)
    lfp.(vn{i}) = cds.lfp.(vn{i});
end
data.lfp = lfp;

units = repmat(struct(),1,length(cds.units));
for i = 1:length(cds.units)
    units(i).id = [cds.units(i).chan, cds.units(i).ID];
    units(i).array = cds.units(i).array;
    units(i).bank = cds.units(i).bank;
    units(i).pin = cds.units(i).pin;
    units(i).label = cds.units(i).label;
    units(i).rowNum = cds.units(i).rowNum;
    units(i).colNum = cds.units(i).colNum;
    units(i).spikes = cds.units(i).spikes;
end
good_idx = [cds.units.ID]~=0 & [cds.units.ID]~=255;
data.units = units(good_idx);

