clear;
clc;

%%
glm_encoding_model_plots;

if 1
    d_vaf = d_vaf./repmat(nanmean(all_m_cv,2),1,size(d_vaf,2));
    ste = mean(nanstd(all_m_cv./repmat(nanmean(all_m_cv,2),1,size(all_m_cv,2)),[],2));
    ylim = [-4,2];
else
    ste = mean(nanstd(all_m_cv,[],2));
    ylim = [-0.15,0.015];
end

a = d_vaf(:,ad_inds);
a = reshape(a,prod(size(a)),1);
b = [ones(length(a),1), reshape(repmat(1:length(ad_inds),size(d_vaf,1),1),length(a),1)];
[fit,~,~,~,s] = regress(a,b);

a = d_vaf(:,wo_inds);
a = reshape(a,prod(size(a)),1);
b = [ones(length(a),1), reshape(repmat(1:length(wo_inds),size(d_vaf,1),1),length(a),1)];
[fit2,~,~,~,s2] = regress(a,b);

%%
close all;

figure;
subplot(1,2,1); hold all;
% plot baseline error for reference
patch([1,201,201,1],[-ste,-ste,ste,ste],[0.7 0.7 0.7])
plot(1:201,nanmean(d_vaf(:,1:201),1),'b+');
plot(1:201,fit(1)+(1:201)*fit(2),'b','LineWidth',3);
axis('tight');
set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',ylim);
title('Force','FontSize',14);

ylabel('Change in pseudo-R^2','FontSize',14);
subplot(1,2,2); hold all;
patch([1,154,154,1],[-ste,-ste,ste,ste],[0.7 0.7 0.7])
plot(1:154,nanmean(d_vaf(:,wo_inds),1),'b+');
plot(1:154,fit2(1)+(1:154)*fit2(2),'b','LineWidth',3);
axis('tight'); set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',ylim);
title('Washout','FontSize',14);

%%
% % % % figure;
% % % subplot(1,2,1); hold all;
% % % plot(1:201,nanmean(d_vaf(:,1:201),1),'r.');
% % % plot(1:201,fit(1)+(1:201)*fit(2),'k','LineWidth',3);
% % % axis('tight');
% % % set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[-3,1]);
% % % title('Force','FontSize',14);
% % % ylabel('Change in pseudo-R^2','FontSize',14);
% % % subplot(1,2,2); hold all;
% % % plot(1:154,nanmean(d_vaf(:,wo_inds),1),'r.');
% % % plot(1:154,fit2(1)+(1:154)*fit2(2),'k','LineWidth',3);
% % % axis('tight'); set(gca,'Box','off','TickDir','out','FontSize',14,'YLim',[-3,1]);
% % % title('Washout','FontSize',14);