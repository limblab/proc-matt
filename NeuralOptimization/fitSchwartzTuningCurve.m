function [tunCurves, sg] = fitSchwartzTuningCurve(filename,holdTime)
% Fit tuning curve to target location during hold period

tau = 0.1; %time in seconds accounting for transmission delay

doCos = true; % convert to cosine fit or leave as sum of sin/cos
doPlots = false;

load(filename);
trialTable = wf_trial_table(out_struct);
neural = out_struct.units;
sg = [neural.id];
sg = sg(1:2:end);

% Using the trial table, compute the hold periods
holdPeriods = [trialTable(:,8)-holdTime, trialTable(:,8)]; % trial end time - hold time

% Get spike count for each channel in each hold period
spikeCounts = zeros(size(holdPeriods,1),length(neural));
for unit = 1:length(neural)
    ts = neural(unit).ts;
    for iTrial = 1:size(holdPeriods,1)
        % how many spikes are in this window?
        spikeCounts(iTrial,unit) = length(ts(ts > holdPeriods(iTrial,1) & ts <= holdPeriods(iTrial,2)));
    end
end

% % Raster code
% figure;
% hold all;
% for unit = 1:length(neural)
%     ts = neural(unit).ts;
%     ts = ts(ts > 100 & ts <= 120);
%     plot([ts'; ts';], [repmat(unit-1,size(ts))'; repmat(unit,size(ts))'],'k');
% end

% Get target data
% Get the target centers
targIDs = trialTable( :, 10);

targetCenters = [(trialTable(:,4)+trialTable(:,2))/2 (trialTable(:,5)+trialTable(:,3))/2];
% map targids to target centers
targIDList = sort(unique(targIDs));
targMap = zeros(length(targIDList),4);
% Get angles and positions to the target centers
for i = 1:length(targIDList)
    target = targetCenters(find(targIDs==targIDList(i),1),:);
    targMap(i,:) = [atan2(target(2),target(1)) targIDList(i) target(1) target(2)];
end


%%% Here is the data we need
fr = spikeCounts./holdTime; % Compute a firing rate
theta = targMap(targIDs,1); % Get angles at each trial's target ID
theta = wrapAngle(theta,0); % make sure it goes from [-pi,pi)

% Now find the tuning curves for each cell
tunCurves = zeros(size(fr,2),3);

st = sin(theta);
ct = cos(theta);
X = [ones(size(theta)) st ct];

for iN = 1:size(fr,2)
    % model is b0+b1*cos(theta)+b2*sin(theta)
    b = regress(fr(:,iN),X);
    if ~doCos
        [~, k] = max(X*b);
        pds(iN) = theta(k);
        if doPlots
            plot(theta,X*b,'b.')
            hold all
            plot(theta,fr(:,iN),'r.')
            % Line to plot PD? May not work
            plot([pds(iN) pds(iN)],[0 20],'k')
            pause;
            close all
        end
    else
        % convert to model b0 + b1*cos(theta+b2)
        b  = [b(1); sqrt(b(2).^2 + b(3).^2); atan2(b(2),b(3))];
        pds(iN) = b(3);
        if doPlots
            temp = b(1) + b(2)*cos(theta-b(3));
            [~,I] = sort(theta);
            plot(theta(I),temp(I),'b','LineWidth',2)
            hold all
            plot(theta,fr(:,iN),'r.')
            plot([pds(iN) pds(iN)],[0 20],'k')
            pause;
            close all
        end
    end
    tunCurves(iN,:) = b;
end