sg = trial_data(1).M1_unit_guide;

channels = 1:96;

cell_counts = zeros(96,1);
for i = channels
    cell_counts(i) = sum(sg(:,1)==i);
end

close all;
figure;
hist(cell_counts,0:max(cell_counts));
ylabel('# Channels','FontSize',14);
xlabel('# Sorted Units','FontSize',14);
title(['Chewie 2016-09-15; Total count = ' num2str(length(sg))],'FontSize',16);
axis tight;
set(gca,'Box','off','TickDir','out','FontSize',14);