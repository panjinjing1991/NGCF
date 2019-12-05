load('/Users/lewlee/Documents/MATLAB/CASE_USA/AI_USA.mat')
load('/Users/lewlee/Documents/MATLAB/CASE_USA/CONUS_mask.mat')
for i = 1:112
    for j = 1:232
        if CONUS_mask_qdeg(i,j) ==0
            AI(113-i,j) = NaN;
        end
    end
end
latt = [25:0.25:52.75]';
lonn = [-125:0.25:-67.25]';

figure;    
%subplot('Position',[0.01 0.51 0.5 0.45]);
subplot('position',[0.15 0.5 0.4 0.4])

m_proj('Lambert Conformal Conic','lat',[25 50],'lon',[-125 -65]);
h=m_pcolor(lonn,latt,flipud(AI));                                 %作图
m_coast('color',[0 0 0],'linewidth',1)                     %画出海岸线
m_gshhs('hb1','linewidth',1,'color','k'); %画出国界线?
hold on;
set(h,'edgecolor','none');
m_grid('box','on','xtick',7,'ytick',6,'fontsize',12,'fontname','Times New Roman','gridcolor','none') 
%m_colmap('jet')
colormap(prism(4))
caxis([0.5 4.5])
colorbar('location','SouthOutside','Ticks',[],'FontSize',15,...
    'XTick',1:1:4,'XtickLabel',{'Dry','Semi-Dry','Semi-Wet','Wet'},...
    'fontname','Times New Roman'); 
title('(a) Climate region','Fontsize',20,'fontname','Times New Roman');


subplot('position',[0.18 0.08 0.35 0.35])
%load('/Users/lewlee/Documents/MATLAB/CASE_USA/Result/spImpact_USA_spatial.mat')
load('/Users/lewlee/Documents/MATLAB/CASE_USA/Result/Case_S_NCA_P_NCA.mat')

load('/Users/lewlee/Documents/MATLAB/CASE_USA/AI_USA.mat')

% wet-day
% read data
Save_pval = Savep;
Select_grd = nan(232,112);
for i = 1:232
    for j = 1:112
        if Save_pval(i,j) < 0.05
            Select_grd(i,j) = 1;
        end
    end
end

Save_S_impact = SavespImpact;

SP_impact_dry = nan(232,112);
SP_impact_wet = nan(232,112);

for i = 1:232
    for j = 1:112
        if Select_grd(i,j) == 1
            SP_impact_dry(i,j) = squeeze(Save_S_impact(1,i,j));
            SP_impact_wet(i,j) = squeeze(Save_S_impact(2,i,j));
        end
    end
end

SP_impact_wet = rot90(SP_impact_wet);
SP_impact_dry = rot90(SP_impact_dry);

y = SP_impact_wet(find(AI == 1));
%y = SP_impact_dry(find(AI == 1));

y(isnan(y)) = [];

ymin=min(y);
ymax=max(y);
% x=linspace(ymin,ymax,100);
% yy=hist(y,x);
% yy=yy/length(y);
% bar(x,yy)

x1=linspace(ymin,ymax,100);
yy1=hist(y,x1);
yy1=yy1/length(y);
f1 = ksdensity(y,x1,'function','pdf');
f1 = f1*(ymax-ymin)./100;
plot(x1,f1,'LineWidth',3)


y = SP_impact_wet(find(AI == 2));
%y = SP_impact_dry(find(AI == 2));

y(isnan(y)) = [];

ymin=min(y);
ymax=max(y);
% x=linspace(ymin,ymax,100);
% yy=hist(y,x);
% yy=yy/length(y);
% bar(x,yy)

x1=linspace(ymin,ymax,100);
yy1=hist(y,x1);
yy1=yy1/length(y);
f1 = ksdensity(y,x1,'function','pdf');
f1 = f1*(ymax-ymin)./100;
hold on
plot(x1,f1,'LineWidth',3)

y = SP_impact_wet(find(AI == 3));
%y = SP_impact_dry(find(AI == 3));

y(isnan(y)) = [];

ymin=min(y);
ymax=max(y);
% x=linspace(ymin,ymax,100);
% yy=hist(y,x);
% yy=yy/length(y);
% bar(x,yy)

x1=linspace(ymin,ymax,100);
yy1=hist(y,x1);
yy1=yy1/length(y);
f1 = ksdensity(y,x1,'function','pdf');
f1 = f1*(ymax-ymin)./100;
hold on
plot(x1,f1,'LineWidth',3)

y = SP_impact_wet(find(AI == 4));
%y = SP_impact_dry(find(AI == 4));

y(isnan(y)) = [];

ymin=min(y);
ymax=max(y);
% x=linspace(ymin,ymax,100);
% yy=hist(y,x);
% yy=yy/length(y);
% bar(x,yy)

x1=linspace(ymin,ymax,100);
yy1=hist(y,x1);
yy1=yy1/length(y);
f1 = ksdensity(y,x1,'function','pdf');
f1 = f1*(ymax-ymin)./100;
hold on
plot(x1,f1,'LineWidth',3)
xlabel('SM-POCC impact','FontSize',16,'fontname','Times New Roman')
ylabel('Probability density','FontSize',16,'fontname','Times New Roman')
legend('Dry','Semi-Dry','Semi-Wet','Wet','fontname','Times New Roman')
title('(b) Probability density curve of SM-P','FontSize',20,'fontname','Times New Roman')
