figure,
subplot(2,3,1),histogram(cross_cor.gR_all,'Normalization','probability')
xlabel('gR (h-1)')
ylabel('fraction of cells')
a=mean(cross_cor.gR_all);
b=std(cross_cor.gR_all);
title(strcat('mean:',num2str(a,'%.2f'),32,'std:',num2str(b,'%.2f'),32 ,'CV:',num2str(b/a,'%.2f')),'FontSize',16);
%32 ´ú±í¿Õ¸ñ

subplot(2,3,2),histogram(cross_cor.fluo1_all,'Normalization','probability')
xlabel('fluo1 (a.u.)')
ylabel('fraction of cells')
a=mean(cross_cor.fluo1_all);
b=std(cross_cor.fluo1_all);
title(strcat('mean:',num2str(a,'%.2f'),32,'std:',num2str(b,'%.2f'),32 ,'CV:',num2str(b/a,'%.2f')),'FontSize',16);

subplot(2,3,3),histogram(cross_cor.p_all,'Normalization','probability')
xlabel('p')
ylabel('fraction of cells')
a=mean(cross_cor.p_all);
b=std(cross_cor.p_all);
title(strcat('mean:',num2str(a,'%.2f'),32,'std:',num2str(b,'%.2f'),32 ,'CV:',num2str(b/a,'%.2f')),'FontSize',16);

subplot(2,3,4),histogram(cross_cor.gRpks_all(2,:),'Normalization','probability')
xlabel('gR pks intensity(h-1)')
ylabel('fraction of cells')
a=mean(cross_cor.gRpks_all(2,:));
b=std(cross_cor.gRpks_all(2,:));
title(strcat('mean:',num2str(a,'%.2f'),32,'std:',num2str(b,'%.2f'),32 ,'CV:',num2str(b/a,'%.2f')),'FontSize',16);

subplot(2,3,5),histogram(cross_cor.fluo1pks_all(2,:),'Normalization','probability')
xlabel('fluo1 pks intensity (a.u.)')
ylabel('fraction of cells')
a=mean(cross_cor.fluo1pks_all(2,:));
b=std(cross_cor.fluo1pks_all(2,:));
title(strcat('mean:',num2str(a,'%.2f'),32,'std:',num2str(b,'%.2f'),32 ,'CV:',num2str(b/a,'%.2f')),'FontSize',16);

subplot(2,3,6),histogram(cross_cor.ppks_all(2,:),'Normalization','probability')
xlabel('p pks intensity (a.u.)')
ylabel('fraction of cells')
a=mean(cross_cor.ppks_all(2,:));
b=std(cross_cor.ppks_all(2,:));
title(strcat('mean:',num2str(a,'%.2f'),32,'std:',num2str(b,'%.2f'),32 ,'CV:',num2str(b/a,'%.2f')),'FontSize',16);

saveas(gca,'hist_analysis','fig');