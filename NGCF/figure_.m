% Figure method of output dataframe in NGCF
% the result dataframe is [impact,p_value,varargout]x[nLon,nLat]
%
%
function fig = figure_()
fig.run_figure = @run_figure;
fig.sign = @sign;
fig.run = @run;
end

function run(result)
figure
[impact_dry,impact_wet] = sign(result);
subplot('Position',[0.01 0.01 0.48 0.45]);
run_figure(impact_dry,'jet','Dry')
subplot('Position',[0.01 0.51 0.48 0.45]);
run_figure(impact_wet,'jet','Wet')
end

function fig = run_figure(X,color,title_)

lat = [25:0.25:52.75]';
lon = [-125:0.25:-67.25]';

m_proj('Lambert Conformal Conic',...
       'lat',[25 50],...
       'lon',[-125 -65]);
h = m_pcolor(lon,lat,X');
m_coast('color',[0 0 0],'linewidth',1);
m_gshhs('hb1','linewidth',1,'color','k'); 
hold on;
set(h,'edgecolor','none');
m_grid('box','on',...
       'xtick',7,...
       'ytick',6,...
       'fontsize',12,...
       'fontname','Times New Roman',...
       'gridcolor','none');
%m_colmap('jet')
colormap(color);
caxis([0.8,1.2]);
colorbar('location','SouthOutside','FontSize',11); 
title(title_,'Fontsize',20);

end

function [impact_dry,impact_wet] = sign(result)
% give p_value shape as (nLon,nLat)

%
impact_dry = squeeze(result(1,:,:));
impact_wet = squeeze(result(2,:,:));
p_value = squeeze(result(3,:,:));
%
[nLon,nLat] = size(p_value);
%
%impact_dry(find(p_value>0.05))= nan;
%impact_wet(find(p_value>0.05))= nan;

end
