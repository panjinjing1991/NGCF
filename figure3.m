GLM = squeeze(Save_cv_R2(:,2,:,:,:));
RF = squeeze(Save_cv_R2(:,1,:,:,:));
RF1 = 1-RF;
GLM1 = 1-GLM;

RF2 = RF1.*aa;
RF3 = squeeze(nanmean(RF2,2));
RF_test = squeeze((RF3(1,:,:)./round(3287*0.2)).^(0.5));
RF_train = squeeze((RF3(2,:,:)./round(3287*0.8)).^(0.5));

GLM2 = GLM1.*aa;
GLM3 = squeeze(nanmean(GLM2,2));
GLM_test = squeeze((GLM3(1,:,:)./round(3287*0.2)).^(0.5));
GLM_train = squeeze((GLM3(2,:,:)./round(3287*0.8)).^(0.5));

latt = [25:0.25:52.75]';
lonn = [-125:0.25:-67.25]';

figure;    
%subplot('Position',[0.01 0.51 0.5 0.45]);
subplot(1,2,1)
m_proj('Lambert Conformal Conic','lat',[25 50],'lon',[-125 -65]);
h=m_pcolor(lonn,latt,GLM_test');                                 %作图
m_coast('color',[0 0 0],'linewidth',1)                     %画出海岸线
m_gshhs('hb1','linewidth',1,'color','k'); %画出国界线?
hold on;
set(h,'edgecolor','none');
m_grid('box','on','xtick',7,'ytick',6,'fontsize',12,'fontname','Times New Roman','gridcolor','none') 
m_colmap('jet')
caxis([0.4,0.6])
colorbar('location','SouthOutside','FontSize',11); 
title('(a) GLM','Fontsize',20);

%subplot('Position',[0.25 0.04 0.5 0.45]);
subplot(1,2,2)
m_proj('Lambert Conformal Conic','lat',[25 50],'lon',[-125 -65]);
h=m_pcolor(lonn,latt,RF_test');                                 %作图
m_coast('color',[0 0 0],'linewidth',1)                     %画出海岸线
m_gshhs('hb1','linewidth',1,'color','k'); %画出国界线?
hold on;
set(h,'edgecolor','none');
m_grid('box','on','xtick',7,'ytick',6,'fontsize',12,'fontname','Times New Roman','gridcolor','none') 
m_colmap('jet')
caxis([0.4,0.6])
colorbar('location','SouthOutside','FontSize',11); 
title('(b) RF','Fontsize',20);

%print(gcf,'-dpdf','/Users/lewlee/Desktop/RMSE.pdf');
print(gcf,'-dtiff','/Users/lewlee/Desktop/RMSE.tiff');


