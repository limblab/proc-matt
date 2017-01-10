function out = fitTuningCurves_Reg(data,params,tuningPeriod,useArray,doPlots)
% notes about inputs
% notes about outputs

if nargin < 6
    doPlots = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load all of the parameters
includeSpeed = params.tuning.includeSpeed;
confLevel = params.tuning.confidenceLevel;
bootNumIters = params.tuning.numberBootIterations;
tuningStatTest = params.tuning.tuningStatTest;

%% Get data
sg = data.(useArray).sg;

[fr,theta,mt,force_rms,vel,force,theta_hand] = getFR(data,params,useArray,tuningPeriod);

% Do bootstrapping with regression
switch lower(tuningStatTest)
    case 'bootstrap'
        statTestParams = {'bootstrap',bootNumIters,confLevel};
    case 'anova'
        statTestParams = {'anova',confLevel};
    case 'none'
        statTestParams = {'none'};
end

% the output will be cell array with element for each block if desired
for iBlock = 1:length(fr)
    disp(['%% Block ' num2str(iBlock) ' of ' num2str(length(fr)) '...']);
    
    if strcmpi(params.paramSetName,'movement_mp') && strcmpi(data.meta.epoch,'ad') && strcmpi(data.meta.perturbation,'ff')
        disp('Using motor plan a la Vaadia...');
        [tcs,cbs,rs,boot_pds,boot_mds,boot_bos] = regressTuningCurves(fr{iBlock},theta_hand{iBlock},statTestParams,'doplots',doPlots);
    else
        if ~includeSpeed
            [tcs,cbs,rs,boot_pds,boot_mds,boot_bos] = regressTuningCurves(fr{iBlock},theta{iBlock},statTestParams,'doplots',doPlots);
        else
            disp('Using Moran/Schwartz speed model...');
            [tcs,cbs,rs,boot_pds,boot_mds,boot_bos] = regressTuningCurves_Moran(fr{iBlock},theta{iBlock},vel{iBlock},statTestParams,'doplots',doPlots);
        end
    end
    

    
    pds = tcs(:,3);
    pd_cis = cbs{3};
    mds = tcs(:,2);
    md_cis = cbs{2};
    bos = tcs(:,1);
    bo_cis = cbs{1};
    
    out(iBlock).pds = [pds pd_cis];
    out(iBlock).mds = [mds md_cis];
    out(iBlock).bos = [bos bo_cis];
    
    out(iBlock).boot_pds = boot_pds;
    out(iBlock).boot_mds = boot_mds;
    out(iBlock).boot_bos = boot_bos;
    out(iBlock).r_squared = rs;
    
    out(iBlock).sg = sg;
    out(iBlock).fr = fr{iBlock};
    out(iBlock).theta = theta{iBlock};
    out(iBlock).theta_hand = theta_hand{iBlock};
    out(iBlock).mt = mt{iBlock};
    out(iBlock).forces = force_rms{iBlock};
    out(iBlock).force_mean = force{iBlock};
    out(iBlock).vels = vel{iBlock};
    out(iBlock).params = params;
    
    out(iBlock).meta = data.meta;
end

