clear;
close all;

baseDir = 'C:\Users\Matt Perich\Desktop\lab\data\';
usePeriod = 'initial';
tuneMethod = 'regression';
targdir = 'target';

binSize = 5;

plotMax = 90;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
useArray = 'M1';

% useDates = {'Chewie','2013-10-28','FF','RT',targdir,'FF'; ...
%             'Chewie','2013-10-29','FF','RT',targdir,'FF'; ...
%             'Chewie','2013-12-09','FF','RT',targdir,'FF'; ...
%             'Chewie','2013-12-10','FF','RT',targdir,'FF'; ...
%             'Chewie','2013-12-12','VR','RT',targdir,'VR'; ...
%             'Chewie','2013-12-13','VR','RT',targdir,'VR'; ...
%             'Chewie','2013-12-17','FF','RT',targdir,'FF'; ... 
%             'Chewie','2013-12-18','FF','RT',targdir,'FF'; ...
%             'Chewie','2013-12-19','VR','CO',targdir,'VR'; ...
%             'Chewie','2013-12-20','VR','CO',targdir,'VR'; ...
%             'Mihili','2014-01-14','VR','RT',targdir,'VR'; ...
%             'Mihili','2014-01-15','VR','RT',targdir,'VR'; ...
%             'Mihili','2014-02-18-VR','VR','CO',targdir,'VR'; ...
%             'Mihili','2014-03-03','VR','CO',targdir,'VR'; ...
%             'Mihili','2014-02-03','FF','CO',targdir,'FF'; ...
%             'Mihili','2014-02-14','FF','RT',targdir,'FF'; ...
%             'Mihili','2014-02-17','FF','CO',targdir,'FF'; ...
%             'Mihili','2014-02-18','FF','CO',targdir,'FF'; ...
%             'Mihili','2014-02-24','FF','RT',targdir,'FF'; ...
%             'Mihili','2014-02-24-VR','VR','RT',targdir,'VR'; ...
%             'Mihili','2014-03-07','FF','CO',targdir,'FF'; ...
%             'Chewie','2013-10-09','VR','RT',targdir,'VR'; ...
%             'Chewie','2013-10-10','VR','RT',targdir,'VR'; ...
%             'Chewie','2013-10-11','VR','RT',targdir,'VR'}; 
useDates = {'Chewie','2013-10-28','FF','RT',targdir,'FF'; ...
            'Chewie','2013-10-29','FF','RT',targdir,'FF'; ...
            'Chewie','2013-12-09','FF','RT',targdir,'FF'; ...
            'Chewie','2013-12-10','FF','RT',targdir,'FF'; ...
            'Chewie','2013-12-12','VR','RT',targdir,'VR'; ...
            'Chewie','2013-12-13','VR','RT',targdir,'VR'; ...
            'Chewie','2013-12-17','FF','RT',targdir,'FF'; ... 
            'Chewie','2013-12-18','FF','RT',targdir,'FF'; ...
            'Mihili','2014-01-14','VR','RT',targdir,'VR'; ...
            'Mihili','2014-01-15','VR','RT',targdir,'VR'; ...
            'Mihili','2014-02-14','FF','RT',targdir,'FF'; ...
            'Mihili','2014-02-24','FF','RT',targdir,'FF'; ...
            'Mihili','2014-02-24-VR','VR','RT',targdir,'VR'; ...
            'Chewie','2013-10-09','VR','RT',targdir,'VR'; ...
            'Chewie','2013-10-10','VR','RT',targdir,'VR'; ...
            'Chewie','2013-10-11','VR','RT',targdir,'VR'}; 
        
        
plotPDShiftComparisonHistograms('dir',baseDir,'dates',useDates,'period',usePeriod,'tunemethod',tuneMethod,'array',useArray,'binsize',binSize,'useblocks',1,'plotmax',plotMax);

useArray = 'PMd';
% useDates = {'Mihili','2014-01-14','VR','RT',targdir,'VR'; ...
%             'Mihili','2014-01-15','VR','RT',targdir,'VR'; ...
%             'Mihili','2014-02-18-VR','VR','CO',targdir,'VR'; ...
%             'Mihili','2014-03-03','VR','CO',targdir,'VR'; ...
%             'Mihili','2014-02-03','FF','CO',targdir,'FF'; ...
%             'Mihili','2014-02-14','FF','RT',targdir,'FF'; ...
%             'Mihili','2014-02-17','FF','CO',targdir,'FF'; ...
%             'Mihili','2014-02-18','FF','CO',targdir,'FF'; ...
%             'Mihili','2014-02-24','FF','RT',targdir,'FF'; ...
%             'Mihili','2014-02-24-VR','VR','RT',targdir,'VR'; ...
%             'Mihili','2014-03-07','FF','CO',targdir,'FF'; ...
%             'MrT','2013-08-19','FF','CO',targdir,'FF'; ...
%             'MrT','2013-08-20','FF','RT',targdir,'FF'; ...
%             'MrT','2013-08-21','FF','CO',targdir,'FF'; ... 
%             'MrT','2013-08-22','FF','RT',targdir,'FF'; ...
%             'MrT','2013-08-23','FF','CO',targdir,'FF'; ...
%             'MrT','2013-08-30','FF','RT',targdir,'FF'; ...
%             'MrT','2013-09-03','VR','CO',targdir,'VR'; ...
%             'MrT','2013-09-04','VR','RT',targdir,'VR'; ...
%             'MrT','2013-09-05','VR','CO',targdir,'VR'; ...
%             'MrT','2013-09-06','VR','RT',targdir,'VR'; ...
%             'MrT','2013-09-09','VR','CO',targdir,'VR'; ...
%             'MrT','2013-09-10','VR','RT',targdir,'VR'}; 
useDates = {'Mihili','2014-01-14','VR','RT',targdir,'VR'; ...
            'Mihili','2014-01-15','VR','RT',targdir,'VR'; ...
            'Mihili','2014-02-14','FF','RT',targdir,'FF'; ...
            'Mihili','2014-02-24','FF','RT',targdir,'FF'; ...
            'Mihili','2014-02-24-VR','VR','RT',targdir,'VR'; ...
            'MrT','2013-08-20','FF','RT',targdir,'FF'; ...
            'MrT','2013-08-22','FF','RT',targdir,'FF'; ...
            'MrT','2013-08-30','FF','RT',targdir,'FF'; ...
            'MrT','2013-09-04','VR','RT',targdir,'VR'; ...
            'MrT','2013-09-06','VR','RT',targdir,'VR'; ...
            'MrT','2013-09-09','VR','CO',targdir,'VR'; ...
            'MrT','2013-09-10','VR','RT',targdir,'VR'}; 
plotPDShiftComparisonHistograms('dir',baseDir,'dates',useDates,'period',usePeriod,'tunemethod',tuneMethod,'array',useArray,'binsize',binSize,'useblocks',3,'plotmax',plotMax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

useArray = 'M1';

% useDates_VR = {'Mihili','2014-01-14','VR','RT',targdir,'VR'; ...
%             'Mihili','2014-01-15','VR','RT',targdir,'VR'; ...
%             'Mihili','2014-02-18-VR','VR','CO',targdir,'VR'; ...
%             'Mihili','2014-03-03','VR','CO',targdir,'VR'; ...
%             'Chewie','2013-10-09','VR','RT',targdir,'VR'; ...
%             'Chewie','2013-10-10','VR','RT',targdir,'VR'; ...
%             'Chewie','2013-10-11','VR','RT',targdir,'VR'; ...
%             'Chewie','2013-12-12','VR','RT',targdir,'VR'; ...
%             'Chewie','2013-12-13','VR','RT',targdir,'VR'; ...
%             'Chewie','2013-12-19','VR','CO',targdir,'VR'; ...
%             'Chewie','2013-12-20','VR','CO',targdir,'VR'}; 
%         
% useDates_FF = {'Mihili','2014-02-03','FF','CO',targdir,'FF'; ...
%             'Mihili','2014-02-14','FF','RT',targdir,'FF'; ...
%             'Mihili','2014-02-17','FF','CO',targdir,'FF'; ...
%             'Mihili','2014-02-18','FF','CO',targdir,'FF'; ...
%             'Chewie','2013-10-28','FF','RT',targdir,'FF'; ...
%             'Chewie','2013-10-29','FF','RT',targdir,'FF'; ...
%             'Chewie','2013-12-09','FF','RT',targdir,'FF'; ...
%             'Chewie','2013-12-10','FF','RT',targdir,'FF'; ...
%             'Chewie','2013-12-17','FF','RT',targdir,'FF'; ... 
%             'Chewie','2013-12-18','FF','RT',targdir,'FF'}; 
useDates_VR = {'Mihili','2014-01-14','VR','RT',targdir,'VR'; ...
            'Mihili','2014-01-15','VR','RT',targdir,'VR'; ...
            'Chewie','2013-10-09','VR','RT',targdir,'VR'; ...
            'Chewie','2013-10-10','VR','RT',targdir,'VR'; ...
            'Chewie','2013-10-11','VR','RT',targdir,'VR'; ...
            'Chewie','2013-12-12','VR','RT',targdir,'VR'; ...
            'Chewie','2013-12-13','VR','RT',targdir,'VR'}; 
        
useDates_FF = {'Mihili','2014-02-14','FF','RT',targdir,'FF'; ...
            'Chewie','2013-10-28','FF','RT',targdir,'FF'; ...
            'Chewie','2013-10-29','FF','RT',targdir,'FF'; ...
            'Chewie','2013-12-09','FF','RT',targdir,'FF'; ...
            'Chewie','2013-12-10','FF','RT',targdir,'FF'; ...
            'Chewie','2013-12-17','FF','RT',targdir,'FF'; ... 
            'Chewie','2013-12-18','FF','RT',targdir,'FF'}; 
        
[VRm_means_bl,VRm_means_ad,VRm_means_wo] = plotPDShiftComparisonHistograms('dir',baseDir,'dates',useDates_VR,'period',usePeriod,'tunemethod',tuneMethod,'array',useArray,'binsize',binSize,'useblocks',1:3);
[FFm_means_bl,FFm_means_ad,FFm_means_wo] = plotPDShiftComparisonHistograms('dir',baseDir,'dates',useDates_FF,'period',usePeriod,'tunemethod',tuneMethod,'array',useArray,'binsize',binSize,'useblocks',1:3);

VRm_means_bl = VRm_means_bl{1};
FFm_means_bl = FFm_means_bl{1};
VRm_means_wo = VRm_means_wo{1};
FFm_means_wo = FFm_means_wo{1};

figure('Position',[200 200 1280 800]);
hold all;
plot(0,0,'bo','LineWidth',2);
plot(0.1,0,'ro','LineWidth',2);
legend({'Curl Field','Visual Rotation'});
plot([0;0],[FFm_means_bl(1,1)+FFm_means_bl(2,1);FFm_means_bl(1,1)-FFm_means_bl(2,1)],'b','LineWidth',2);
plot([0.1;0.1],[VRm_means_bl(1,1)+VRm_means_bl(2,1);VRm_means_bl(1,1)-VRm_means_bl(2,1)],'r','LineWidth',2);

for iBlock = 1:length(VRm_means_ad)
    tempVR = VRm_means_ad{iBlock};
    tempFF = FFm_means_ad{iBlock};
    
    plot(iBlock,tempFF(1,:),'bo','LineWidth',2);
    plot(iBlock+0.1,tempVR(1,:),'ro','LineWidth',3);
    plot([iBlock;iBlock],[tempFF(1,:)+tempFF(2,:);tempFF(1,:)-tempFF(2,:)],'b','LineWidth',2);
    plot([iBlock+0.1;iBlock+0.1],[tempVR(1,:)+tempVR(2,:);tempVR(1,:)-tempVR(2,:)],'r','LineWidth',2);
end


plot(4:5,FFm_means_wo(1,:),'bo','LineWidth',2);
plot(4.1:5.1,VRm_means_wo(1,:),'ro','LineWidth',3);
plot([4:5;4:5],[FFm_means_wo(1,:)+FFm_means_wo(2,:);FFm_means_wo(1,:)-FFm_means_wo(2,:)],'b','LineWidth',2);
plot([4.1:5.1;4.1:5.1],[VRm_means_wo(1,:)+VRm_means_wo(2,:);VRm_means_wo(1,:)-VRm_means_wo(2,:)],'r','LineWidth',2);

plot([-1 7],[0 0],'LineWidth',1,'Color',[0.6 0.6 0.6]);
plot([0.5 0.5],[-plotMax plotMax],'k--','LineWidth',1);
plot([3.5 3.5],[-plotMax plotMax],'k--','LineWidth',1);


axis([-0.3 5.3 -plotMax plotMax]);


set(gca,'XTick',[0 1 2 3 4 5],'XTickLabel',{'Base','Early','Mid','Late','Early Wash', 'Late Wash'},'FontSize',14);

title([useArray '-' targdir],'FontSize',16);
ylabel('Change in PD (deg)','FontSize',16);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

useArray = 'PMd';

useDates_VR = {'Mihili','2014-01-14','VR','RT',targdir,'VR'; ...
            'Mihili','2014-01-15','VR','RT',targdir,'VR'; ...
            'Mihili','2014-02-18-VR','VR','CO',targdir,'VR'; ...
            'Mihili','2014-03-03','VR','CO',targdir,'VR'; ...
            'MrT','2013-09-03','VR','CO',targdir,'VR'; ...
            'MrT','2013-09-04','VR','RT',targdir,'VR'; ...
            'MrT','2013-09-05','VR','CO',targdir,'VR'; ...
            'MrT','2013-09-06','VR','RT',targdir,'VR'; ...
            'MrT','2013-09-09','VR','CO',targdir,'VR'; ...
            'MrT','2013-09-10','VR','RT',targdir,'VR'}; 
        
useDates_FF = {'Mihili','2014-02-03','FF','CO',targdir,'FF'; ...
            'Mihili','2014-02-14','FF','RT',targdir,'FF'; ...
            'Mihili','2014-02-17','FF','CO',targdir,'FF'; ...
            'Mihili','2014-02-18','FF','CO',targdir,'FF'; ...
            'MrT','2013-08-19','FF','CO',targdir,'FF'; ...
            'MrT','2013-08-20','FF','RT',targdir,'FF'; ...
            'MrT','2013-08-21','FF','CO',targdir,'FF'; ... 
            'MrT','2013-08-22','FF','RT',targdir,'FF'; ...
            'MrT','2013-08-23','FF','CO',targdir,'FF'; ...
            'MrT','2013-08-30','FF','RT',targdir,'FF'}; 

[VRm_means_bl,VRm_means_ad,VRm_means_wo] = plotPDShiftComparisonHistograms('dir',baseDir,'dates',useDates_VR,'period',usePeriod,'tunemethod',tuneMethod,'array',useArray,'binsize',binSize,'useblocks',1:3);
[FFm_means_bl,FFm_means_ad,FFm_means_wo] = plotPDShiftComparisonHistograms('dir',baseDir,'dates',useDates_FF,'period',usePeriod,'tunemethod',tuneMethod,'array',useArray,'binsize',binSize,'useblocks',1:3);

VRm_means_bl = VRm_means_bl{1};
FFm_means_bl = FFm_means_bl{1};
VRm_means_wo = VRm_means_wo{1};
FFm_means_wo = FFm_means_wo{1};

figure('Position',[200 200 1280 800]);
hold all;
plot(0,0,'bo','LineWidth',2);
plot(0.1,0,'ro','LineWidth',2);
legend({'Curl Field','Visual Rotation'});
plot([0;0],[FFm_means_bl(1,1)+FFm_means_bl(2,1);FFm_means_bl(1,1)-FFm_means_bl(2,1)],'b','LineWidth',2);
plot([0.1;0.1],[VRm_means_bl(1,1)+VRm_means_bl(2,1);VRm_means_bl(1,1)-VRm_means_bl(2,1)],'r','LineWidth',2);

for iBlock = 1:length(VRm_means_ad)
    tempVR = VRm_means_ad{iBlock};
    tempFF = FFm_means_ad{iBlock};
    
    plot(iBlock,tempFF(1,:),'bo','LineWidth',2);
    plot(iBlock+0.1,tempVR(1,:),'ro','LineWidth',3);
    plot([iBlock;iBlock],[tempFF(1,:)+tempFF(2,:);tempFF(1,:)-tempFF(2,:)],'b','LineWidth',2);
    plot([iBlock+0.1;iBlock+0.1],[tempVR(1,:)+tempVR(2,:);tempVR(1,:)-tempVR(2,:)],'r','LineWidth',2);
end


plot(4:5,FFm_means_wo(1,:),'bo','LineWidth',2);
plot(4.1:5.1,VRm_means_wo(1,:),'ro','LineWidth',3);
plot([4:5;4:5],[FFm_means_wo(1,:)+FFm_means_wo(2,:);FFm_means_wo(1,:)-FFm_means_wo(2,:)],'b','LineWidth',2);
plot([4.1:5.1;4.1:5.1],[VRm_means_wo(1,:)+VRm_means_wo(2,:);VRm_means_wo(1,:)-VRm_means_wo(2,:)],'r','LineWidth',2);

plot([-1 7],[0 0],'LineWidth',1,'Color',[0.6 0.6 0.6]);
plot([0.5 0.5],[-plotMax plotMax],'k--','LineWidth',1);
plot([3.5 3.5],[-plotMax plotMax],'k--','LineWidth',1);


axis([-0.3 5.3 -plotMax plotMax]);


set(gca,'XTick',[0 1 2 3 4 5],'XTickLabel',{'Base','Early','Mid','Late','Early Wash', 'Late Wash'},'FontSize',14);

title([useArray '-' targdir],'FontSize',16);
ylabel('Change in PD (deg)','FontSize',16);


% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NOW DO ADAPTING CELL COMPARISON
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
useArray = 'M1';
useDates = {'Mihili','2014-01-14','VR','RT',targdir,'VR'; ...
            'Mihili','2014-01-15','VR','RT',targdir,'VR'; ...
            'Mihili','2014-02-18-VR','VR','CO',targdir,'VR'; ...
            'Mihili','2014-03-03','VR','CO',targdir,'VR'; ...
            'Mihili','2014-02-03','FF','CO',targdir,'FF'; ...
            'Mihili','2014-02-14','FF','RT',targdir,'FF'; ...
            'Mihili','2014-02-17','FF','CO',targdir,'FF'; ...
            'Mihili','2014-02-18','FF','CO',targdir,'FF'; ...
            'Chewie','2013-10-09','VR','RT',targdir,'VR'; ...
            'Chewie','2013-10-10','VR','RT',targdir,'VR'; ...
            'Chewie','2013-10-11','VR','RT',targdir,'VR'; ...
            'Chewie','2013-10-28','FF','RT',targdir,'FF'; ...
            'Chewie','2013-10-29','FF','RT',targdir,'FF'; ...
            'Chewie','2013-12-09','FF','RT',targdir,'FF'; ...
            'Chewie','2013-12-10','FF','RT',targdir,'FF'; ...
            'Chewie','2013-12-12','VR','RT',targdir,'VR'; ...
            'Chewie','2013-12-13','VR','RT',targdir,'VR'; ...
            'Chewie','2013-12-17','FF','RT',targdir,'FF'; ... 
            'Chewie','2013-12-18','FF','RT',targdir,'FF'; ...
            'Chewie','2013-12-19','VR','CO',targdir,'VR'; ...
            'Chewie','2013-12-20','VR','CO',targdir,'VR'}; 
[ap,np] = plotAdaptingCellComparison('dir',baseDir,'dates',useDates,'period',usePeriod,'array',useArray,'tunemethod',tuneMethod,'useblocks',3);
title('M1');

useArray = 'PMd';
useDates = {'Mihili','2014-01-14','VR','RT',targdir,'VR'; ...
            'Mihili','2014-01-15','VR','RT',targdir,'VR'; ...
            'Mihili','2014-02-18-VR','VR','CO',targdir,'VR'; ...
            'Mihili','2014-03-03','VR','CO',targdir,'VR'; ...
            'Mihili','2014-02-03','FF','CO',targdir,'FF'; ...
            'Mihili','2014-02-14','FF','RT',targdir,'FF'; ...
            'Mihili','2014-02-17','FF','CO',targdir,'FF'; ...
            'Mihili','2014-02-18','FF','CO',targdir,'FF'; ...
            'MrT','2013-08-19','FF','CO',targdir,'FF'; ...
            'MrT','2013-08-20','FF','RT',targdir,'FF'; ...
            'MrT','2013-08-21','FF','CO',targdir,'FF'; ... 
            'MrT','2013-08-22','FF','RT',targdir,'FF'; ...
            'MrT','2013-08-23','FF','CO',targdir,'FF'; ...
            'MrT','2013-08-30','FF','RT',targdir,'FF'; ...
            'MrT','2013-09-03','VR','CO',targdir,'VR'; ...
            'MrT','2013-09-04','VR','RT',targdir,'VR'; ...
            'MrT','2013-09-05','VR','CO',targdir,'VR'; ...
            'MrT','2013-09-06','VR','RT',targdir,'VR'; ...
            'MrT','2013-09-09','VR','CO',targdir,'VR'; ...
            'MrT','2013-09-10','VR','RT',targdir,'VR'}; 
[ap,np] = plotAdaptingCellComparison('dir',baseDir,'dates',useDates,'period',usePeriod,'array',useArray,'tunemethod',tuneMethod,'useblocks',3);
title('PMd');

% % xInds = [1, 1.1, 1.22, 1.32];
% xInds = [1, 1.1, 1.22];
% figure('Position',[200 200 800 600]);
% hold all;
% 
% colors = {'b','b','b','b'};
% for i = 1:length(xInds)
%     temp = nan(length(xInds),2);
%     temp(i,1) = ap(i);
%     temp(i,2) = np(i);
%     h = bar(xInds,temp,0.8,'stacked');
%     
%     set(get(h(1),'Children'),'FaceColor',colors{i},'EdgeColor',colors{i},'LineWidth',6);
%     set(get(h(2),'Children'),'FaceColor','w','EdgeColor',colors{i},'LineWidth',6);
% end
% 
% axis('tight');
% set(gca,'XTick',xInds,'XTickLabel',{'VR 1','VR 2','VR 3','VR 4'},'FontSize',14);
% title([useArray '-' targdir],'FontSize',16);
% ylabel('Percent of Adapting Cells','FontSize',14);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % make plots showing percentage of adapting cells in early vs late
% [ap,np] = plotAdaptingCellComparison('dir',baseDir,'dates',useDates,'period',usePeriod,'array',useArray,'tunemethod',tuneMethod,'useblocks',3);
% 
% % xInds = [1, 1.1, 1.22, 1.32];
% xInds = [1, 1.1, 1.22];
% figure('Position',[200 200 800 600]);
% hold all;
% 
% colors = {'b','b','b','b'};
% for i = 1:length(xInds)
%     temp = nan(length(xInds),2);
%     temp(i,1) = ap(i);
%     temp(i,2) = np(i);
%     h = bar(xInds,temp,0.8,'stacked');
%     
%     set(get(h(1),'Children'),'FaceColor',colors{i},'EdgeColor',colors{i},'LineWidth',6);
%     set(get(h(2),'Children'),'FaceColor','w','EdgeColor',colors{i},'LineWidth',6);
% end
% 
% axis('tight');
% set(gca,'XTick',xInds,'XTickLabel',{'VR 1','VR 2','VR 3','VR 4'},'FontSize',14);
% title([useArray '-' targdir],'FontSize',16);
% ylabel('Percent of Adapting Cells','FontSize',14);
% 
