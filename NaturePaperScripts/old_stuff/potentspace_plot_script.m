function potentspace_plot_script(td)
close all;
figure;

epochs = {'BL','AD'};
arrow_size = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_idx = 2:50;
marker_loc = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputArgs = struct( ...
    'which_one', 'PMdM1_null', ...
    'axis_range',[-1, 1], ...
    'sp_num',1, ...
    't_idx',t_idx, ...
    'marker_loc',marker_loc, ...
    'arrow_size',arrow_size);
do_the_plot(td,epochs,inputArgs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputArgs = struct( ...
    'which_one', 'PMdM1_potent', ...
    'axis_range',[-1, 1], ...
    'sp_num',2, ...
    't_idx',t_idx, ...
    'marker_loc',marker_loc, ...
    'arrow_size',arrow_size);
do_the_plot(td,epochs,inputArgs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputArgs = struct( ...
    'which_one', 'M1_pca', ...
    'axis_range',[-1, 1.4], ...
    'sp_num',3, ...
    't_idx',t_idx, ...
    'marker_loc',marker_loc, ...
    'arrow_size',arrow_size+2);
do_the_plot(td,epochs,inputArgs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_idx = 65:100;
marker_loc = 70;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputArgs = struct( ...
    'which_one', 'PMdM1_null', ...
    'axis_range',[-1, 1], ...
    'sp_num',1, ...
    't_idx',t_idx, ...
    'marker_loc',marker_loc, ...
    'arrow_size',arrow_size);
do_the_plot(td,epochs,inputArgs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputArgs = struct( ...
    'which_one', 'PMdM1_potent', ...
    'axis_range',[-1, 1], ...
    'sp_num',2, ...
    't_idx',t_idx, ...
    'marker_loc',marker_loc, ...
    'arrow_size',4);
do_the_plot(td,epochs,inputArgs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputArgs = struct( ...
    'which_one', 'M1_pca', ...
    'axis_range',[-1, 1.4], ...
    'sp_num',3, ...
    't_idx',t_idx, ...
    'marker_loc',marker_loc, ...
    'arrow_size',arrow_size+2);
do_the_plot(td,epochs,inputArgs);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function do_the_plot(td,epochs,inputArgs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
which_one = '';
axis_range = [];
sp_num = 0;
t_idx = [];
marker_loc = 0;
arrow_size = 0;
assignParams(who,inputArgs); %get parameters

sp_cols = 3;
sp_rows = 1;

epoch_sym = {'-','--'};
epoch_widths = [1,2];
axis_min = axis_range(1);
axis_max = axis_range(2);

u = unique([td.target_direction]);
c = phasemap(length(u));

subplot(sp_rows,sp_cols,sp_num);
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
title(which_one,'FontSize',14);
axis('square');

set(gca,'Box','off','TickDir','out','FontSize',14);
xlabel('PC1','FontSize',14);
ylabel('PC2','FontSize',14);
zlabel('PC3','FontSize',14);
grid on

end