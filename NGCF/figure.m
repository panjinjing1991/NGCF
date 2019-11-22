function fig = figure()

% figure

%figure('position',[50 50 800 600])

axesm('mercator','maplatlimit',[25.125 52.875],'maplonlimit',[-124.875 -67.125], 'frame', 'off', ...
      'MeridianLabel', 'off', 'ParallelLabel', 'off', ...
      'plinelocation', 25.125:0.25:52.875, 'mlinelocation', -124.875:0.25:-67.125, ...
      'grid', 'off', 'PLabelMeridian', 'west', 'MLabelParallel','south','labelformat','none');

%Z = rot90(Select_grd,3);
subplot(2,1,1)

Z = fliplr(rot90(SP_impact_dry,3));
R = georasterref('RasterSize', size(Z),'Latlim', [25.125 52.875], 'Lonlim', [-124.875 -67.125]);
h = grid2image(Z,R);
set(h,'alphadata',~isnan(Z));
set(gca, 'visible', 'on');
geoshow([UsShape.Y],[UsShape.X],'Color','black')
colormap(jet(8))
caxis([0.8,1.2])
%title('GLM Spatial Dry','color','b','FontSize',15)
title('lag4 dry','color','b','FontSize',15)

colorbar

subplot(2,1,2)

Z = fliplr(rot90(SP_impact_wet,3));
R = georasterref('RasterSize', size(Z),'Latlim', [25.125 52.875], 'Lonlim', [-124.875 -67.125]);
h = grid2image(Z,R);
set(h,'alphadata',~isnan(Z));
set(gca, 'visible', 'on');
geoshow([UsShape.Y],[UsShape.X],'Color','black')
colormap(jet(8))
caxis([0.8,1.2])
%title('GLM Spatial Wet','color','b','FontSize',15)
title('lag4 wet','color','b','FontSize',15)

colorbar



%subplot('Position',[0.01 0.51 0.5 0.45]);
subplot('Position',[0.01 0.01 0.48 0.45]);
m_proj('Lambert Conformal Conic','lat',[25 50],'lon',[-125 -65]);
h=m_pcolor(lonn,latt,SP_impact_dry');                                 %作图
m_coast('color',[0 0 0],'linewidth',1)                     %画出海岸线
m_gshhs('hb1','linewidth',1,'color','k'); %画出国界线?
hold on;
set(h,'edgecolor','none');
m_grid('box','on','xtick',7,'ytick',6,'fontsize',12,'fontname','Times New Roman','gridcolor','none') 
%m_colmap('jet')
colormap('jet')
caxis([0.8,1.2])
colorbar('location','SouthOutside','FontSize',11); 
title('1x1(Dry)','Fontsize',20);

%subplot('Position',[0.25 0.04 0.5 0.45]);
subplot('Position',[0.51 0.01 0.48 0.45]);
m_proj('Lambert Conformal Conic','lat',[25 50],'lon',[-125 -65]);
h=m_pcolor(lonn,latt,SP_impact_wet');                                 %作图
m_coast('color',[0 0 0],'linewidth',1)                     %画出海岸线
m_gshhs('hb1','linewidth',1,'color','k'); %画出国界线?
hold on;
set(h,'edgecolor','none');
m_grid('box','on','xtick',7,'ytick',6,'fontsize',12,'fontname','Times New Roman','gridcolor','none') 
%m_colmap('jet')
colormap('jet')
caxis([0.8,1.2])
colorbar('location','SouthOutside','FontSize',11); 
title('1x1(Wet)','Fontsize',20);

print(gcf,'-djpeg','-r600','/Users/lewlee/Desktop/Spatial.jpg');


end
