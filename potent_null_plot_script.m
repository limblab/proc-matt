%%
close all;
figure;

epochs = {'BL','AD'};
epoch_sym = {'-','--'};
epoch_widths = [1,1];

u = unique([td.target_direction]);
% u = u(2:2:end);
c = phasemap(length(u));

t_idx = 2:50;
marker_loc = 20;

arrow_size = 4;

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % which_one = 'PMd_pca';
% % % axis_min = -1.4;
% % % axis_max = 1.4;
% % % % c = [0 0 1];
% % % p1=subplot(131);
% % % 
% % % for i = 1:length(u)
% % %     for e = 1:length(epochs)
% % %     idx =  getTDidx(td,'target_direction',u(i),'epoch',epochs{e});
% % %     for j = 1:length(idx)
% % %         trial = idx(j);
% % %         hold all;
% % %         temp = td(trial).(which_one);
% % %         plot3(temp(t_idx,1),temp(t_idx,2),temp(t_idx,3),epoch_sym{e},'LineWidth',epoch_widths(e),'Color',c(i,:))
% % %         plot3(temp(marker_loc,1),temp(marker_loc,2),temp(marker_loc,3),'o','LineWidth',1,'Color',c(i,:))
% % %         arrowMMC3(temp(t_idx(end-2),1:3),temp(t_idx(end-1),1:3),temp(t_idx(end),1:3),arrow_size, ...
% % %             [axis_min,axis_max,axis_min,axis_max,axis_min,axis_max],c(i,:),c(i,:));
% % %     end
% % %     end
% % % end
% % % title('PMd','FontSize',14);
% % % axis('square');
% % % set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[axis_min,axis_max],'YLim',[axis_min,axis_max],'ZLim',[axis_min,axis_max]);
% % % xlabel('PC1','FontSize',14);
% % % ylabel('PC2','FontSize',14);
% % % zlabel('PC3','FontSize',14);
% % % grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1=subplot(121);
axis_min = -0.4;
axis_max = 0.4;

for i = 1:length(u)
    for e = 1:length(epochs)
    idx =  getTDidx(td,'target_direction',u(i),'epoch',epochs{e});
    for j = 1:length(idx)
        trial = idx(j);
        hold all;
        temp1 = td(trial).potent;
        temp2 = td(trial).null;
        plot3(temp1(t_idx,1),temp1(t_idx,2),temp2(t_idx,1),epoch_sym{e},'LineWidth',epoch_widths(e),'Color',c(i,:))
        plot3(temp1(marker_loc,1),temp1(marker_loc,2),temp2(marker_loc,1),'o','LineWidth',1,'Color',c(i,:))
        arrowMMC3([temp1(t_idx(end-2),1:2),temp2(t_idx(end-2),1)], ...
            [temp1(t_idx(end-1),1:2),temp2(t_idx(end-1),1)], ...
            [temp1(t_idx(end),1:2),temp2(t_idx(end),1)], ...
            arrow_size, ...
            [axis_min,axis_max,axis_min,axis_max,axis_min,axis_max], ...
            c(i,:),c(i,:));
    end
    end
end
title('Potent/Null','FontSize',14);
axis('square');
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[axis_min,axis_max],'YLim',[axis_min,axis_max],'ZLim',[axis_min,axis_max]);
xlabel('Potent 1','FontSize',14);
ylabel('Potent 2','FontSize',14);
zlabel('Null 1','FontSize',14);
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1=subplot(122);
arrow_size = 6;
which_one = 'M1_pca';
axis_min = -1.4;
axis_max = 1.4;

for i = 1:length(u)
    for e = 1:length(epochs)
    idx =  getTDidx(td,'target_direction',u(i),'epoch',epochs{e});
    for j = 1:length(idx)
        trial = idx(j);
        hold all;
        temp = td(trial).(which_one);
        plot3(temp(t_idx,1),temp(t_idx,2),temp(t_idx,3),epoch_sym{e},'LineWidth',epoch_widths(e),'Color',c(i,:))
        plot3(temp(marker_loc,1),temp(marker_loc,2),temp(marker_loc,3),'o','LineWidth',1,'Color',c(i,:))
        arrowMMC3(temp(t_idx(end-2),1:3),temp(t_idx(end-1),1:3),temp(t_idx(end),1:3),arrow_size, ...
            [axis_min,axis_max,axis_min,axis_max,axis_min,axis_max],c(i,:),c(i,:));
    end
    end
end
title('M1 PCA','FontSize',14);
axis('square');
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[axis_min,axis_max],'YLim',[axis_min,axis_max],'ZLim',[axis_min,axis_max]);
xlabel('PC1','FontSize',14);
ylabel('PC2','FontSize',14);
zlabel('PC3','FontSize',14);
grid on



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_idx = 65:100;
marker_loc = 70;
arrow_size = 4;
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % which_one = 'PMd_pca';
% % % axis_min = -1.4;
% % % axis_max = 1.4;
% % % % c = [0 0 1];
% % % n2=subplot(131);
% % % 
% % % for i = 1:length(u)
% % %     for e = 1:length(epochs)
% % %     idx =  getTDidx(td,'target_direction',u(i),'epoch',epochs{e});
% % %     for j = 1:length(idx)
% % %         trial = idx(j);
% % %         hold all;
% % %         temp = td(trial).(which_one);
% % %         plot3(temp(t_idx,1),temp(t_idx,2),temp(t_idx,3),epoch_sym{e},'LineWidth',epoch_widths(e),'Color',c(i,:))
% % %         plot3(temp(marker_loc,1),temp(marker_loc,2),temp(marker_loc,3),'o','LineWidth',1,'Color',c(i,:))
% % %         arrowMMC3(temp(t_idx(end-2),1:3),temp(t_idx(end-1),1:3),temp(t_idx(end),1:3),arrow_size, ...
% % %             [axis_min,axis_max,axis_min,axis_max,axis_min,axis_max],c(i,:),c(i,:));
% % %     end
% % %     end
% % % end
% % % title('PMd','FontSize',14);
% % % axis('square');
% % % set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[axis_min,axis_max],'YLim',[axis_min,axis_max],'ZLim',[axis_min,axis_max]);
% % % xlabel('PC1','FontSize',14);
% % % ylabel('PC2','FontSize',14);
% % % zlabel('PC3','FontSize',14);
% % % grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p2=subplot(121);
axis_min = -0.4;
axis_max = 0.4;

for i = 1:length(u)
    for e = 1:length(epochs)
    idx =  getTDidx(td,'target_direction',u(i),'epoch',epochs{e});
    for j = 1:length(idx)
        trial = idx(j);
        hold all;
        temp1 = td(trial).potent;
        temp2 = td(trial).null;
        plot3(temp1(t_idx,1),temp1(t_idx,2),temp2(t_idx,1),epoch_sym{e},'LineWidth',epoch_widths(e),'Color',c(i,:))
        plot3(temp1(marker_loc,1),temp1(marker_loc,2),temp2(marker_loc,1),'o','LineWidth',1,'Color',c(i,:))
        arrowMMC3([temp1(t_idx(end-2),1:2),temp2(t_idx(end-2),1)], ...
            [temp1(t_idx(end-1),1:2),temp2(t_idx(end-1),1)], ...
            [temp1(t_idx(end),1:2),temp2(t_idx(end),1)], ...
            arrow_size, ...
            [axis_min,axis_max,axis_min,axis_max,axis_min,axis_max], ...
            c(i,:),c(i,:));
    end
    end
end
title('Potent/Null','FontSize',14);
axis('square');
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[axis_min,axis_max],'YLim',[axis_min,axis_max],'ZLim',[axis_min,axis_max]);
xlabel('Potent 1','FontSize',14);
ylabel('Potent 2','FontSize',14);
zlabel('Null 1','FontSize',14);
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m2=subplot(122);
arrow_size = 6;
which_one = 'M1_pca';
axis_min = -1;
axis_max = 1.4;

for i = 1:length(u)
    for e = 1:length(epochs)
    idx =  getTDidx(td,'target_direction',u(i),'epoch',epochs{e});
    for j = 1:length(idx)
        trial = idx(j);
        hold all;
        temp = td(trial).(which_one);
        plot3(temp(t_idx,1),temp(t_idx,2),temp(t_idx,3),epoch_sym{e},'LineWidth',epoch_widths(e),'Color',c(i,:))
        plot3(temp(marker_loc,1),temp(marker_loc,2),temp(marker_loc,3),'o','LineWidth',1,'Color',c(i,:))
        arrowMMC3(temp(t_idx(end-2),1:3),temp(t_idx(end-1),1:3),temp(t_idx(end),1:3),arrow_size, ...
            [axis_min,axis_max,axis_min,axis_max,axis_min,axis_max],c(i,:),c(i,:));

    end
    end
end
title('M1 PCA','FontSize',14);
axis('square');
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[axis_min,axis_max],'YLim',[axis_min,axis_max],'ZLim',[axis_min,axis_max]);
xlabel('PC1','FontSize',14);
ylabel('PC2','FontSize',14);
zlabel('PC3','FontSize',14);
grid on
