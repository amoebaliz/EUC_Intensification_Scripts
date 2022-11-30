close all; clear all; clc

load mean_fields_500.mat
load redblue.mat


I_lat = find(latx>=-0.5&latx<=0.5);
depth = 1:40;
axislims_1 = [150 270 -1 1];
axislims_2 = [150 270 9 34];

figure
set(gcf, 'Renderer', 'Painters', 'Position', [-1469 -86 970 971])

%% frich

hsub = subplot(2,2,1);
%set(hsub, 'Position',[.57 .35 .34 .14])

[~,h]=contourf(lon,depth,(10^7)*(squeeze(mean(frich5(:,I_lat,:),2))'),100);
set(h,'EdgeColor','none')
hold on

contour(lon,1:40,squeeze(mean(u5(:,I_lat,:),2))',[0.5 0.5],'k','LineWidth',1);
axis(axislims_2)

% Y axis
set(gca,'YDir','reverse')

set(gca,'TickDir','out','Ytick',[10 20 30],'fontsize', 14,'fontweight', 'bold')

% X axis

set(gca,'TickDir','out','Xtick',[],'XTickLabel',{[]})

% C axis
cmax = .2;
hcb = colorbar; colormap(redblue); caxis([-1*cmax cmax]);

set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
    'fontweight', 'normal')

%set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
%    'fontweight', 'normal','Position',[0.95 0.358 0.016 0.102])

%set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
%    'fontname','arial','fontsize',12,'Position',  [1.6 1.3*cmax 1])

set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
    'fontname','arial','fontsize',12)

% plot border
set(gca,'box','off')
isholdonque = ishold; 
hold on
ax = axis;
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
plot(ax(1:2),(ax(3))*[1,1],'k','linewidth',1)
if isholdonque == 0
hold off
end 

% plot headings
text(150, 6,'\textsf{\textbf{f)\ }}{\boldmath$\frac{\partial}{\partial{x}}\bigl[A_{H}\frac{\partial{u}}{\partial{x}}\bigr] + \frac{\partial}{\partial{y}}\bigl[A_{H}\frac{\partial{u}}{\partial{y}}\bigr]$}',...
    'interpreter','latex','FontSize',14,'fontweight', 'bold', 'fontname','arial');

%% fricv

hsub = subplot(2,2,3);
%set(hsub, 'Position',[.57 .16 .34 .14])

[~,h]=contourf(lon,depth,(10^7)*(squeeze(mean(fricv5(:,I_lat,:),2))'),100);
set(h,'EdgeColor','none')
hold on

contour(lon,1:40,squeeze(mean(u5(:,I_lat,:),2))',[0.5 0.5],'k','LineWidth',1);
axis(axislims_2)

% Y axis
set(gca,'YDir','reverse')

set(gca,'TickDir','out','Ytick',[10 20 30],'fontsize', 14,'fontweight', 'bold')

% X axis
set(gca,'TickDir','out','Xtick',[150 210 270],...
   'XTickLabel',{'150 E', '150 W','90 W'},...
   'fontsize', 14,'fontweight','bold')

xlabel('Longitude','fontweight', 'bold', 'Position', [210 42 1]);

% C axis
cmax = 5;
hcb = colorbar; colormap(redblue); caxis([-1*cmax cmax]);


set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
    'fontweight', 'normal')

%set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
%    'fontweight', 'normal','Position',[0.95 0.168 0.016 0.102])

%set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
%    'fontname','arial','fontsize',12,'Position',  [1.6 1.3*cmax 1])

set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
    'fontname','arial','fontsize',12)



% plot border
set(gca,'box','off')
isholdonque = ishold; 
hold on
ax = axis;
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
plot(ax(1:2),(ax(3))*[1,1],'k','linewidth',1)
if isholdonque == 0
hold off
end 

% plot headings
text(150, 6,'\textsf{\textbf{h)\ }}{\boldmath$\frac{\partial}{\partial{z}}\bigl[A_{V}\frac{\partial{u}}{\partial{z}}\bigr]$}',...
    'interpreter','latex','FontSize',14,'fontweight', 'bold', 'fontname','arial');

%%
load mean_fields_8000.mat

%% frich

hsub = subplot(2,2,2);
% set(hsub, 'Position',[.57 .35 .34 .14])

[~,h]=contourf(lon,depth,(10^7)*(squeeze(mean(frich5(:,I_lat,:),2))'),100);
set(h,'EdgeColor','none')
hold on

contour(lon,1:40,squeeze(mean(u5(:,I_lat,:),2))',[0.5 0.5],'k','LineWidth',1);
axis(axislims_2)

% Y axis
set(gca,'YDir','reverse')

set(gca,'TickDir','out','Ytick',[10 20 30],'fontsize', 14,'fontweight', 'bold')

% X axis

set(gca,'TickDir','out','Xtick',[],'XTickLabel',{[]})

% C axis
cmax = .2;
hcb = colorbar; colormap(redblue); caxis([-1*cmax cmax]);
set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
    'fontweight', 'normal')

% set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
%    'fontweight', 'normal','Position',[0.95 0.358 0.016 0.102])

% set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
%    'fontname','arial','fontsize',12,'Position',  [1.6 1.3*cmax 1])

set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
    'fontname','arial','fontsize',12)

% plot border
set(gca,'box','off')
isholdonque = ishold; 
hold on
ax = axis;
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
plot(ax(1:2),(ax(3))*[1,1],'k','linewidth',1)
if isholdonque == 0
hold off
end 

% plot headings
text(150, 6,'\textsf{\textbf{f)\ }}{\boldmath$\frac{\partial}{\partial{x}}\bigl[A_{H}\frac{\partial{u}}{\partial{x}}\bigr] + \frac{\partial}{\partial{y}}\bigl[A_{H}\frac{\partial{u}}{\partial{y}}\bigr]$}',...
    'interpreter','latex','FontSize',14,'fontweight', 'bold', 'fontname','arial');

%% fricv

hsub = subplot(2,2,4);
%set(hsub, 'Position',[.57 .16 .34 .14])

[~,h]=contourf(lon,depth,(10^7)*(squeeze(mean(fricv5(:,I_lat,:),2))'),100);
set(h,'EdgeColor','none')
hold on

contour(lon,1:40,squeeze(mean(u5(:,I_lat,:),2))',[0.5 0.5],'k','LineWidth',1);
axis(axislims_2)

% Y axis
set(gca,'YDir','reverse')

set(gca,'TickDir','out','Ytick',[10 20 30],'fontsize', 14,'fontweight', 'bold')

% X axis
set(gca,'TickDir','out','Xtick',[150 210 270],...
   'XTickLabel',{'150 E', '150 W','90 W'},...
   'fontsize', 14,'fontweight','bold')

xlabel('Longitude','fontweight', 'bold', 'Position', [210 42 1]);

% C axis
cmax = 5;
hcb = colorbar; colormap(redblue); caxis([-1*cmax cmax]);
set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
    'fontweight', 'normal')

%set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
%    'fontweight', 'normal','Position',[0.95 0.168 0.016 0.102])

%set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
%    'fontname','arial','fontsize',12,'Position',  [1.6 1.3*cmax 1])

set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
    'fontname','arial','fontsize',12)



% plot border
set(gca,'box','off')
isholdonque = ishold; 
hold on
ax = axis;
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
plot(ax(1:2),(ax(3))*[1,1],'k','linewidth',1)
if isholdonque == 0
hold off
end 

% plot headings
text(150, 6,'\textsf{\textbf{h)\ }}{\boldmath$\frac{\partial}{\partial{z}}\bigl[A_{V}\frac{\partial{u}}{\partial{z}}\bigr]$}',...
    'interpreter','latex','FontSize',14,'fontweight', 'bold', 'fontname','arial');