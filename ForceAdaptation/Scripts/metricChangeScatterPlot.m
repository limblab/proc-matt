% This is for a general "population PD change scatter plot". Horizontal and
% vertical axis can be specified independently to be:
%   - Any epoch value or epoch difference (eg BL->AD, BL->WO, etc)
%   - Any metric (PD, MD, BO, FR)
%   - Any coordinate (e.g. movement vs target)
%
%   This is all mediated with the toCompare struct. This code assumes that
%   only the first two entries in that struct array are relevant.

do_dir_anova = true;

classColors = {'k','b','r','m','g'};
classNames = {'Kinematic','Dynamic','Memory I','Memory II','Other'};

%% Hardcode define axis information for plots
metricInfo.PD.min = -180;
metricInfo.PD.max = 180;
metricInfo.PD.binSize = 5;
metricInfo.PD.label = 'Preferred Direction (Deg) ';

metricInfo.MD.min = 0;
metricInfo.MD.max = 50;
metricInfo.MD.binSize = 2;
metricInfo.MD.label = 'Modulation Depth (Hz) ';

metricInfo.BO.min = 0;
metricInfo.BO.max = 50;
metricInfo.BO.binSize = 2;
metricInfo.BO.label = 'Offset (Hz) ';

metricInfo.FR.min = 0;
metricInfo.FR.max = 50;
metricInfo.FR.binSize = 2;
metricInfo.FR.label = 'Mean Firing Rate (Hz) ';

metricInfo.dPD.min = -180;
metricInfo.dPD.max = 180;
metricInfo.dPD.binSize = 5;
metricInfo.dPD.label = 'PD Change (Deg) ';

if sComp.doPercent
    metricInfo.dMD.min = -1.1;
    metricInfo.dMD.max = 1.1;
    metricInfo.dMD.binSize = 0.05;
    metricInfo.dMD.label = 'MD Change';
    
    metricInfo.dBO.min = -1;
    metricInfo.dBO.max = 1;
    metricInfo.dBO.binSize = 0.05;
    metricInfo.dBO.label = 'BO Change';
    
    metricInfo.dFR.min = -1;
    metricInfo.dFR.max = 1;
    metricInfo.dFR.binSize = 0.05;
    metricInfo.dFR.label = 'FR Change';
else
    metricInfo.dMD.min = -30;
    metricInfo.dMD.max = 30;
    metricInfo.dMD.binSize = 2;
    metricInfo.dMD.label = 'MD Change (Hz) ';
    
    metricInfo.dBO.min = -30;
    metricInfo.dBO.max = 30;
    metricInfo.dBO.binSize = 2;
    metricInfo.dBO.label = 'BO Change (Hz) ';
    
    metricInfo.dFR.min = -30;
    metricInfo.dFR.max = 30;
    metricInfo.dFR.binSize = 2;
    metricInfo.dFR.label = 'FR Change (Hz) ';
end

%%
axisNames = {'x','y'};


%%
% If we want to separate by waveform width...
if doWidthSeparation
    count = 0;
    clear allWFWidths;
    for iFile = 1:size(doFiles,1)
        % load baseline data to get width of all spike waveforms
        data = loadResults(root_dir,doFiles(iFile,:),'data',[],'BL');
        
        units = data.(useArray).units;
        for u = 1:length(units)
            count = count + 1;
            wf = mean(units(u).wf,2);
            idx = find(abs(wf) > std(wf));
            allWFWidths(count) = idx(end) - idx(1);
        end
    end
    % now, set the threshold for narrow and wide APs
    wfThresh = median(allWFWidths);
end

%% Get the master spike guide list
for iAxis = 1:length(axisNames)
    paramSetName = sComp.params{iAxis};
    useArray = sComp.arrays{iAxis};
    tuneMethod = sComp.methods{iAxis};
    tuneWindow = sComp.windows{iAxis};
    for iFile = 1:size(doFiles,1)
        % load tuning and class info
        [t,c] = loadResults(root_dir,doFiles(iFile,:),'tuning',{'tuning','classes'},useArray,paramSetName,tuneMethod,tuneWindow);

        classifierBlocks = c(whichBlock).params.classes.classifierBlocks;
        
        if doWidthSeparation
            % load baseline data to get waveforms
            data = loadResults(root_dir,doFiles(iFile,:),'data',[],'BL');
            
            units = data.(useArray).units;
            fileWidths = zeros(length(units),1);
            wfTypes = zeros(length(units),1);
            for u = 1:length(units)
                wf = mean(units(u).wf,2);
                idx = find(abs(wf) > std(wf));
                switch doWidthSeparation
                    case 1
                        wfTypes(u) = (idx(end) - idx(1)) <= wfThresh;
                    case 2
                        wfTypes(u) = (idx(end) - idx(1)) > wfThresh;
                    case 3
                        wfTypes(u) = 1;
                end
                fileWidths(u) = idx(end) - idx(1);
            end
        else
            wfTypes = ones(size(c(whichBlock).istuned,1),1);
            fileWidths = ones(size(c(whichBlock).istuned,1),1);
        end
        tunedCells = c(whichBlock).sg(all(c(whichBlock).istuned(:,whichTuned),2) & ~all(c(whichBlock).istuned,2) & wfTypes,:);
        
        sg_bl = t(classifierBlocks(1)).sg;
        sg_ad = t(classifierBlocks(2)).sg;
        sg_wo = t(classifierBlocks(3)).sg;
        
        [~,idx_bl] = intersect(sg_bl, tunedCells,'rows');
        [~,idx_ad] = intersect(sg_ad, tunedCells,'rows');
        [~,idx_wo] = intersect(sg_wo, tunedCells,'rows');
    end
end

%% Get the classification for each day for tuned cells
for iAxis = 1:length(axisNames)
    paramSetName = sComp.params{iAxis};
    useArray = sComp.arrays{iAxis};
    tuneMethod = sComp.methods{iAxis};
    tuneWindow = sComp.windows{iAxis};
    
    cellClasses = cell(size(doFiles,1),1);
    cellWidths = cell(size(doFiles,1),1);
    cellPDs = cell(size(doFiles,1),1);
    cellMDs = cell(size(doFiles,1),1);
    cellBOs = cell(size(doFiles,1),1);
    cellFRs = cell(size(doFiles,1),1);
    cellR2 = cell(size(doFiles,1),1);
    cellSG = cell(size(doFiles,1),1);
    count = 0;
    
    all_p = [];
    
    pertDir = zeros(1,size(doFiles,1));
    for iFile = 1:size(doFiles,1)
        % load tuning and class info
        [t,c] = loadResults(root_dir,doFiles(iFile,:),'tuning',{'tuning','classes'},useArray,paramSetName,tuneMethod,tuneWindow);
        
        
        if do_dir_anova
            disp('Doing f-test on direction...');
            fr = t(1).fr; theta = t(1).theta;
            p = zeros(1,size(fr,2));
            for k = 1:size(fr,2)
                p(k) = anovan(fr(:,k),theta,'display','off');
            end
        end
        
        
        
        
        % get direction of perturbation to flip the clockwise ones to align
        if flipClockwisePerts
            pdDir = '';
            % gotta hack it
            dataPath = fullfile(root_dir,doFiles{iFile,1},pdDir,doFiles{iFile,2});
            expParamFile = fullfile(dataPath,[doFiles{iFile,2} '_experiment_parameters.dat']);
            t(1).params.exp = parseExpParams(expParamFile);
            switch lower(t(1).params.exp.angle_dir)
                case 'ccw'
                pertDir(iFile) = 1;
                case 'cw'
                pertDir(iFile) = -1;
            end
        else
            pertDir(iFile) = 1;
        end
        
        classifierBlocks = c(whichBlock).params.classes.classifierBlocks;
        
        if doWidthSeparation
            % load baseline data to get waveforms
            data = loadResults(root_dir,doFiles(iFile,:),'data',[],'BL');
            
            units = data.(useArray).units;
            fileWidths = zeros(length(units),1);
            wfTypes = zeros(length(units),1);
            for u = 1:length(units)
                wf = mean(units(u).wf,2);
                idx = find(abs(wf) > std(wf));
                switch doWidthSeparation
                    case 1
                        wfTypes(u) = (idx(end) - idx(1)) <= wfThresh;
                    case 2
                        wfTypes(u) = (idx(end) - idx(1)) > wfThresh;
                    case 3
                        wfTypes(u) = 1;
                end
                fileWidths(u) = idx(end) - idx(1);
            end
        else
            wfTypes = ones(size(c(whichBlock).istuned,1),1);
            fileWidths = ones(size(c(whichBlock).istuned,1),1);
        end
        
        % first column is PD, second column is MD
        tuned_idx = all(c(whichBlock).istuned(:,whichTuned),2) & wfTypes;
        
        if sComp.highSNR
            data = loadResults(root_dir,doFiles(iFile,:),'data','M1','BL');
            cell_snr = zeros(size(data.sg,1),1);
            for unit = 1:size(data.sg,1)
                wf = data.units(unit).wf;
                %istuned(unit,1) = max(rms(wf'))/mean(std(double(wf(1:5,:)'))) >= waveSNR;
                sig = max(mean(wf'))-min(mean(wf'));
                cell_snr(unit) = sig / (2*mean(std(double(wf'))));
            end
            tuned_idx = tuned_idx & cell_snr > 6;
        end
        cellClasses{iFile} = c(whichBlock).classes(tuned_idx,1);
        tunedCells = c(whichBlock).sg(tuned_idx,:);
        
        sg_bl = t(classifierBlocks(1)).sg;
        sg_ad = t(classifierBlocks(2)).sg;
        sg_wo = t(classifierBlocks(3)).sg;
        
        [~,idx_bl] = intersect(sg_bl, tunedCells,'rows');
        [~,idx_ad] = intersect(sg_ad, tunedCells,'rows');
        [~,idx_wo] = intersect(sg_wo, tunedCells,'rows');
        
        pd_bl = t(classifierBlocks(1)).pds(idx_bl,1);
        pd_ad = t(classifierBlocks(2)).pds(idx_ad,1);
        pd_wo = t(classifierBlocks(3)).pds(idx_wo,1);
        
        cellPDs{iFile} = {pd_bl, pd_ad, pd_wo};
        
        md_bl = t(classifierBlocks(1)).mds(idx_bl,1);
        md_ad = t(classifierBlocks(2)).mds(idx_ad,1);
        md_wo = t(classifierBlocks(3)).mds(idx_wo,1);
        
        cellMDs{iFile} = {md_bl, md_ad, md_wo};
        
        bo_bl = t(classifierBlocks(1)).bos(idx_bl,1);
        bo_ad = t(classifierBlocks(2)).bos(idx_ad,1);
        bo_wo = t(classifierBlocks(3)).bos(idx_wo,1);
        
        cellBOs{iFile} = {bo_bl, bo_ad, bo_wo};
        
        fr_bl = mean(t(classifierBlocks(1)).fr(:,idx_bl),1)';
        fr_ad = mean(t(classifierBlocks(2)).fr(:,idx_ad),1)';
        fr_wo = mean(t(classifierBlocks(3)).fr(:,idx_wo),1)';
        
        cellFRs{iFile} = {fr_bl, fr_ad, fr_wo};
        
        r2_bl = mean(t(classifierBlocks(1)).r_squared,2);
        r2_ad = mean(t(classifierBlocks(2)).r_squared,2);
        r2_wo = mean(t(classifierBlocks(3)).r_squared,2);
        
        cellWidths{iFile} = fileWidths(all(c(whichBlock).istuned(:,whichTuned),2) & wfTypes);
        cellR2{iFile} = {r2_bl(idx_bl), r2_ad(idx_ad), r2_wo(idx_wo)};
        
        cellSG{iFile} = sg_bl(idx_bl,:);
        
        if do_dir_anova
            all_p = [all_p, p(all(c(whichBlock).istuned(:,[1,3]),2))];
        end
    end
    
    %
    % find all dPDs
    pd_bl = []; pd_ad = []; pd_wo = [];
    pd_bl_ad = []; pd_bl_wo = []; pd_ad_wo = [];
    md_bl = []; md_ad = []; md_wo = [];
    md_bl_ad = []; md_bl_wo = []; md_ad_wo = [];
    bo_bl = []; bo_ad = []; bo_wo = [];
    bo_bl_ad = []; bo_bl_wo = []; bo_ad_wo = [];
    fr_bl = []; fr_ad = []; fr_wo = [];
    fr_bl_ad = []; fr_bl_wo = []; fr_ad_wo = [];
    classes = [];
    r2s = [];
    widths = [];
    all_pds = [];
    all_sg = [];
    for iFile = 1:size(doFiles,1)
        pds = cellPDs{iFile};
        mds = cellMDs{iFile};
        bos = cellBOs{iFile};
        frs = cellFRs{iFile};
        
        sg = cellSG{iFile};
        
        all_pds = [all_pds; pds{1}];
        all_sg = [all_sg; sg];
        
        pd_bl = [pd_bl; pds{1}.*(180/pi)];
        pd_ad = [pd_ad; pds{2}.*(180/pi)];
        pd_wo = [pd_wo; pds{3}.*(180/pi)];
        
        md_bl = [md_bl; mds{1}];
        md_ad = [md_ad; mds{2}];
        md_wo = [md_wo; mds{3}];
        
        bo_bl = [bo_bl; bos{1}];
        bo_ad = [bo_ad; bos{2}];
        bo_wo = [bo_wo; bos{3}];
        
        fr_bl = [fr_bl; frs{1}];
        fr_ad = [fr_ad; frs{2}];
        fr_wo = [fr_wo; frs{3}];
        
        if sComp.doPercent
            md_mean = mds{1};
            bo_mean = bos{1};
            fr_mean = frs{1};
        else
            md_mean = ones(size(mds{1},1),1);
            bo_mean = ones(size(bos{1},1),1);
            fr_mean = ones(size(frs{1},1),1);
        end
        
        % pertDir is perturbation direction... 1=counter-clockwise, -1=clockwise... so this flips clockwise perturbations to be counter-clockwise
        pd_bl_ad = [pd_bl_ad; pertDir(iFile)*angleDiff(pds{1},pds{2},true,true).*(180/pi)];
        pd_bl_wo = [pd_bl_wo; pertDir(iFile)*angleDiff(pds{1},pds{3},true,true).*(180/pi)];
        pd_ad_wo = [pd_ad_wo; pertDir(iFile)*angleDiff(pds{2},pds{3},true,true).*(180/pi)];
        
        md_bl_ad = [md_bl_ad; (mds{2}-mds{1})./md_mean];
        md_bl_wo = [md_bl_wo; (mds{3}-mds{1})./md_mean];
        md_ad_wo = [md_ad_wo; (mds{3}-mds{2})./md_mean];
        
        bo_bl_ad = [bo_bl_ad; (bos{2}-bos{1})./bo_mean];
        bo_bl_wo = [bo_bl_wo; (bos{3}-bos{1})./bo_mean];
        bo_ad_wo = [bo_ad_wo; (bos{3}-bos{2})./bo_mean];
        
        fr_bl_ad = [fr_bl_ad; (frs{2}-frs{1})./fr_mean];
        fr_bl_wo = [fr_bl_wo; (frs{3}-frs{1})./fr_mean];
        fr_ad_wo = [fr_ad_wo; (frs{3}-frs{2})./fr_mean];
        
        c = cellClasses{iFile};
        w = cellWidths{iFile};
        
        if sComp.reassignOthers
            pd_wo = pds{3};
            pd_ad = pds{2};
            pd_bl = pds{1};
            idx = find(c==5);
            
            for i = 1:length(idx)
                % calculate the memory cell index
                bl_wo = pertDir(iFile)*angleDiff( pd_wo(idx(i),1),pd_bl(idx(i),1),true,true );
                bl_ad = pertDir(iFile)*angleDiff( pd_ad(idx(i),1), pd_bl(idx(i),1),true,true );
                ad_wo = pertDir(iFile)*angleDiff(pd_wo(idx(i),1),pd_ad(idx(i),1),true,true);
                
                mem_ind = abs(bl_wo) / min( abs(bl_ad) , abs(ad_wo) );
                
                % we also want both BL->WO and BL->AD to be same direction for memory
                %   otherwise it's just dynamic
                if mem_ind > 1
                    if sign(bl_wo)==sign(bl_ad)
                        c(idx(i)) = 3;
                    else
                        c(idx(i)) = 2;
                    end
                elseif mem_ind < 1
                    c(idx(i)) = 2;
                else
                    disp('Hey! This one is exactly one.');
                    c(idx(i)) = 3;
                end
            end
        end
        
        classes = [classes; c];
        widths = [widths; w];
        r2 = cellR2{iFile};
        r2s = [r2s; r2{1}];
    end
    
    % take absolute value of differences if requested
    if sComp.doAbs
        pd_bl_ad = abs(pd_bl_ad); pd_bl_wo = abs(pd_bl_wo); pd_ad_wo = abs(pd_ad_wo);
        md_bl_ad = abs(md_bl_ad); md_bl_wo = abs(md_bl_wo); md_ad_wo = abs(md_ad_wo);
        bo_bl_ad = abs(bo_bl_ad); bo_bl_wo = abs(bo_bl_wo); bo_ad_wo = abs(bo_ad_wo);
        fr_bl_ad = abs(fr_bl_ad); fr_bl_wo = abs(fr_bl_wo); fr_ad_wo = abs(fr_ad_wo);
    end
    
    % see if class correlates with R2
    % plot(r2s,classes,'o');
    
    % decide what data to plot on two axes based on toCompare variable
    if length(sComp.epochs{iAxis}) > 2
        temp = 'd';
    else
        temp = '';
    end
    eval([axisNames{iAxis} '_data = ' lower(sComp.metrics{iAxis}) '_' lower(sComp.epochs{iAxis}) ';']);
    eval([axisNames{iAxis} 'min = metricInfo.([ temp sComp.metrics{iAxis} ]).min;']);
    eval([axisNames{iAxis} 'max = metricInfo.([ temp sComp.metrics{iAxis} ]).max;']);
    eval([axisNames{iAxis} 'BinSize = metricInfo.([ temp sComp.metrics{iAxis} ]).binSize;']);
    eval([axisNames{iAxis} '_lab = [metricInfo.([ temp sComp.metrics{iAxis} ]).label sComp.epochs{iAxis} ];']);
    
end

if sComp.reassignOthers
    classColors = classColors(1:end-1);
    classNames = classNames(1:end-1);
end

figure('Position',[300 50 950 950]);
subplot1(2,2,'Gap',[0 0],'Min',[0.08 0.08],'Max',[0.98 0.98]);

% make a legend here
subplot1(2);

set(gca,'XTick',[],'YTick',[],'XLim',[-0.1 4],'YLim',[-0.1 length(classColors)+1]);
V = axis(gca);
patch([V(1) V(2) V(2) V(1)],[V(3) V(3) V(4) V(4)],[0.8 0.8 0.8]);

for i=1:length(classColors)
    plot(1,i,'d','Color',classColors{i},'LineWidth',2);
    text(2,i,classNames{i},'FontSize',14);
end

box off;

xbins = xmin:xBinSize:xmax;
ybins = ymin:yBinSize:ymax;

% plot histogram of changes BL->AD
subplot1(1);
hist(x_data,xbins);

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','facealpha',1,'edgealpha',1);

set(gca,'YTick',[],'XTick',[],'XLim',[xmin,xmax]);
hold all;
V = axis(gca);

m = mean(x_data);
plot([m m],V(3:4),'k','LineWidth',2);

box off;

% add mean to histogram

% now plot BL->WO histogram
subplot1(4);
hist(y_data,ybins);

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','facealpha',1,'edgealpha',1);

hold all;
V = axis(gca);
m = mean(y_data);
plot([m m],V(3:4),'k','LineWidth',2);

set(gca,'XTick',[],'YTick',[],'XLim',[ymin,ymax]);
set(gca,'CameraUpVector',[1,0,0],'Xdir','reverse');

box off;

% plot scatter of PD changes
subplot1(3);
hold all;
for i = 1:length(x_data)
    plot(x_data(i),y_data(i),'d','LineWidth',2,'Color',classColors{classes(i)});
end

set(gca,'XLim',[xmin,xmax],'YLim',[ymin,ymax],'FontSize',14,'TickDir','out')

V = axis(gca);
plot([V(1) V(2)],[0 0],'k--');
plot([0 0],[V(3) V(4)],'k--');
% plot([V(1) V(2)],[V(3) V(4)],'k-');

plot([mean(x_data) mean(x_data)],V(3:4),'k','LineWidth',2);
plot(V(1:2),[mean(y_data) mean(y_data)],'k','LineWidth',2);
xlabel([x_lab 'BL->AD'],'FontSize',14);
ylabel(y_lab,'FontSize',14);

box off;




