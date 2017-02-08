close all;
clear;
clc;

% chewie array: bundle exits along sulcus laterally
% wire bundle exits from "top", so y=10
% thus, based on normal view, Chewie's 1 on horizontal is against sulcus
% thus, flip Mihili along the horizontal direction
useMonkeys = {'Chewie'};
runAnalysisScripts;
% [all_sg,I] = unique(all_sg,'rows');
% classes = classes(I);

cmp_file = 'Z:\lab_folder\Animal-Miscellany\Chewie 8I2\Blackrock implant surgery 6-14-10\1025-0394.cmp';
% cmp_file = 'Z:\lab_folder\Animal-Miscellany\Mihili 12A3\M1 SN  6250-000989.cmp';

cmp = read_cmp(cmp_file);
cmp = cell2mat(cmp(:,1:3));
% 
mem_counts_chewie = zeros(10);
kin_counts_chewie = zeros(10);
dyn_counts_chewie = zeros(10);
for i = 1:size(cmp,1)
    mem_counts_chewie(cmp(i,1)+1,cmp(i,2)+1) = sum(all_sg(:,1) == cmp(i,3) & (classes == 3 | classes == 4 | classes == 5));
    kin_counts_chewie(cmp(i,1)+1,cmp(i,2)+1) = sum(all_sg(:,1) == cmp(i,3) & (classes == 1));
    dyn_counts_chewie(cmp(i,1)+1,cmp(i,2)+1) = sum(all_sg(:,1) == cmp(i,3) & (classes == 2));
end

% mihili array: bundle exits along sulcus towards midline
useMonkeys = {'Mihili'};
runAnalysisScripts;
% [all_sg,I] = unique(all_sg,'rows');
% classes = classes(I);

cmp_file = 'Z:\lab_folder\Animal-Miscellany\Mihili 12A3\M1 SN  6250-000989.cmp';

cmp = read_cmp(cmp_file);
cmp = cell2mat(cmp(:,1:3));
% 
mem_counts_mihili = zeros(10);
kin_counts_mihili = zeros(10);
dyn_counts_mihili = zeros(10);
for i = 1:size(cmp,1)
    mem_counts_mihili(cmp(i,1)+1,cmp(i,2)+1) = sum(all_sg(:,1) == cmp(i,3) & (classes == 3 | classes == 4 | classes == 5));
    kin_counts_mihili(cmp(i,1)+1,cmp(i,2)+1) = sum(all_sg(:,1) == cmp(i,3) & (classes == 1));
    dyn_counts_mihili(cmp(i,1)+1,cmp(i,2)+1) = sum(all_sg(:,1) == cmp(i,3) & (classes == 2));
end

% flip mihili
mem_counts_mihili = flipud(fliplr(mem_counts_mihili));
kin_counts_mihili = flipud(fliplr(kin_counts_mihili));
dyn_counts_mihili = flipud(fliplr(dyn_counts_mihili));

mem_counts = mem_counts_chewie | mem_counts_mihili;
kin_counts = kin_counts_chewie | kin_counts_mihili;
dyn_counts = dyn_counts_chewie | dyn_counts_mihili;

close all;
figure;

subplot(1,3,1);
imagesc(1:10,1:10,dyn_counts);
axis square;
set(gca,'Box','off','TickDir','out','FontSize',14,'YDir','normal');
title('Dynamic cells','FontSize',14);
% colorbar('SouthOutside');

subplot(1,3,2);
imagesc(1:10,1:10,kin_counts);
axis square;
set(gca,'Box','off','TickDir','out','FontSize',14,'YDir','normal');
title('Kinematic cells','FontSize',14);
% colorbar('SouthOutside');

subplot(1,3,3);
imagesc(1:10,1:10,mem_counts);
axis square;
set(gca,'Box','off','TickDir','out','FontSize',14,'YDir','normal');
title('Memory cells','FontSize',14);
% colorbar('SouthOutside');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now look at dpd over array


