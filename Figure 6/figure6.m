close all; clear all; clc

load mean_fields_E.mat
load redblue.mat


I_lat = find(latx>=-0.5&latx<=0.5);
depth = 1:40;
axislims_1 = [150 270 -1 1];
axislims_2 = [150 270 9 34];

figure
set(gcf, 'Renderer', 'Painters', 'Position', [-1469 -86 970 971])

% Taux and SSH

taux = squeeze(mean(taux_mean(:,I_lat),2));
ssh = squeeze(mean(ssh_mean(:,I_lat),2));

cmax = 1;
yrange = (linspace(-1*cmax,cmax,length(redblue)))';
%
subh = subplot(4,2,1);

set(subh, 'Position', [.09 .73 .34 .14])


%
for j = 1:length(taux)
   
    I_1 = find(yrange <= taux(j),1,'last');
    I_2 = find(yrange >= taux(j),1,'first');
    
    diff_val = taux(j) - yrange(I_1);
    diff_val2 = yrange(I_2)- yrange(I_1);
    frac_val = (diff_val/diff_val2);
    
    taux_colorval = ((1-frac_val)*redblue(I_1,:))+((frac_val)*redblue(I_2,:));
 
    I_1 = find(yrange <= ssh(j),1,'last');
    I_2 = find(yrange >= ssh(j),1,'first');
    
    diff_val = ssh(j) - yrange(I_1);
    diff_val2 = yrange(I_2)- yrange(I_1);
    frac_val = (diff_val/diff_val2);
    
    ssh_colorval = (1-frac_val)*redblue(I_1,:)+(frac_val)*redblue(I_2,:);
    
    plot(lon(j), taux(j),'Color',taux_colorval,'LineWidth',2)
    hold on
    plot(lon(j), ssh(j),'Color',ssh_colorval,'LineWidth',2)
    axis(axislims_1)
    
end

% Y axis
set(gca,'TickDir','out','YTick',[-1 0 1],'fontsize', 14)

ylabel('dynes cm^{-2} & m', 'FontSize',14,'fontweight', 'bold',...
    'Position',[132 .03 1])

% X axis

set(gca,'TickDir','out','Xtick',[],'XTickLabel',{[]})

% Colorbar
hcb = colorbar; colormap(redblue); caxis([-1*cmax cmax])

set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
    'fontweight', 'normal','Position',[0.47 0.738 0.016 0.102])

set(get(hcb,'title'),'string','(dynes cm^{-2} & m)',...
    'fontname','arial','fontsize',12,'Position',  [2.5 1.3*cmax 1])

% plot border
set(gca,'box','off')
isholdonque = ishold; 
hold on
ax = axis;
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
plot(ax(1:2),(ax(4))*[1,1],'k','linewidth',1)
if isholdonque == 0
hold off
end 

% plot headings

text(150, 1.25,'\textsf{\textbf{a)\ }}{\boldmath$\tau^{x}$} \textbf{\&\ SSH}',...
    'interpreter','latex','FontSize',14,'fontweight', 'bold', 'fontname','arial');


%% -1*ududx

hsub = subplot(4,2,3);
set(hsub, 'Position',[.09 .54 .34 .14])

[~,h]=contourf(lon,depth,10^7*(squeeze(mean(-1*ududx5(:,I_lat,:),2))'),100);
set(h,'EdgeColor','none')
hold on

contour(lon,1:40,squeeze(mean(u5(:,I_lat,:),2))',[0.5 0.5],'k','LineWidth',1);
axis(axislims_2)

% Y axis
set(gca,'YDir','reverse')

set(gca,'TickDir','out','Ytick',[10 20 30],'fontsize', 14,'fontweight', 'bold')
% X axis

set(gca,'TickDir','out','Xtick',[])

% C axis
cmax = 5;
hcb = colorbar; colormap(redblue); caxis([-1*cmax cmax]);

set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
    'fontweight', 'normal','Position',[0.47 0.548 0.016 0.102])

set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
    'fontname','arial','fontsize',12,'Position',  [1.6 1.3*cmax 1])


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
text(150, 6,'\textsf{\textbf{c)\ -\ }}{\boldmath$u\frac{\partial{u}}{\partial{x}}$}','interpreter','latex',...
    'FontSize',14,'fontweight', 'bold', 'fontname','arial');


%% -1*vdudy

hsub = subplot(4,2,5);
set(hsub, 'Position',[.09 .35 .34 .14])

[~,h]=contourf(lon,depth,(10^7)*(squeeze(mean(-1*vdudy5(:,I_lat,:),2))'),100);
set(h,'EdgeColor','none')
hold on

contour(lon,1:40,squeeze(mean(u5(:,I_lat,:),2))',[0.5 0.5],'k','LineWidth',1);
axis(axislims_2)

% Y axis
set(gca,'YDir','reverse')

set(gca,'TickDir','out','Ytick',[10 20 30],'fontsize', 14,'fontweight', 'bold')

ylabel('Isopycnal layer','fontweight', 'bold', 'Position', [130 20 1],'fontsize', 16,'fontweight', 'bold')

% X axis

set(gca,'TickDir','out','Xtick',[],'XTickLabel',{[]})

% C axis
cmax = 5;
hcb = colorbar; colormap(redblue); caxis([-1*cmax cmax]);
set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
    'fontweight', 'normal','Position',[0.47 0.358 0.016 0.102])

set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
    'fontname','arial','fontsize',12,'Position',  [1.6 1.3*cmax 1])

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
text(150, 6,'\textsf{\textbf{e)\ -\ }}{\boldmath$v\frac{\partial{u}}{\partial{y}}$}','interpreter','latex',...
    'FontSize',14,'fontweight', 'bold', 'fontname','arial');

%% -1*wdudz

hsub = subplot(4,2,7);
set(hsub, 'Position',[.09 .16 .34 .14])

[~,h]=contourf(lon,depth,(10^7)*(squeeze(mean(-1*wdudz5(:,I_lat,:),2))'),100);
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
    'fontweight', 'normal','Position',[0.47 0.168 0.016 0.102])

set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
    'fontname','arial','fontsize',12,'Position',  [1.6 1.3*cmax 1])


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
text(150, 6,'\textsf{\textbf{g)\ -\ }}{\boldmath$w\frac{\partial{u}}{\partial{z}}$}','interpreter','latex',...
    'FontSize',14,'fontweight', 'bold', 'fontname','arial');


%% u
subh = subplot(4,2,2);

set(subh, 'Position', [.57 .73 .34 .14])

[~,h]=contourf(lon,depth,squeeze(mean(u5(:,I_lat,:),2))',100);
set(h,'EdgeColor','none')
hold on

contour(lon,1:40,squeeze(mean(u5(:,I_lat,:),2))',[0.5 0.5],'k','LineWidth',1);
axis(axislims_2)


% Y axis
set(gca,'YDir','reverse')

set(gca,'TickDir','out','Ytick',[10 20 30],'fontsize', 14,'fontweight', 'bold')

% X axis

set(gca,'TickDir','out','Xtick',[],'XTickLabel',{[]})

% Colorbar
cmax = 1;
hcb = colorbar; colormap(redblue); caxis([-1*cmax cmax])

set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
    'fontweight', 'normal','Position',[0.95 0.738 0.016 0.102])

set(get(hcb,'title'),'string','(m s^{-1})',...
    'fontname','arial','fontsize',12,'Position',  [1.6 1.3*cmax 1])

% plot border
set(gca,'box','off')
isholdonque = ishold; 
hold on
ax = axis;
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
plot(ax(1:2),((ax(3)))*[1,1],'k','linewidth',1)
if isholdonque == 0
hold off
end 

% plot headings
text(150, 6,'\textsf{\textbf{b)\ }}{\boldmath$u$}','interpreter','latex',...
    'FontSize',14,'fontweight', 'bold');

%% zpgf
hsub = subplot(4,2,4);
set(hsub, 'Position',[.57 .54 .34 .14])

[~,h]=contourf(lon,depth,10^7*(squeeze(mean(zpgf5(:,I_lat,:),2))'),100);
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
cmax = 5;
hcb = colorbar; colormap(redblue); caxis([-1*cmax cmax]);

set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
    'fontweight', 'normal','Position',[0.95 0.548 0.016 0.102])

set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
    'fontname','arial','fontsize',12,'Position',  [1.6 1.3*cmax 1])


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
text(150, 6,'\textsf{\textbf{d)\ -\ }}{\boldmath$\frac{1}{\rho}\frac{\partial{P}}{\partial{x}}$}','interpreter','latex',...
    'FontSize',14,'fontweight', 'bold', 'fontname','arial');

%% frich

hsub = subplot(4,2,6);
set(hsub, 'Position',[.57 .35 .34 .14])

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
cmax = 5;
hcb = colorbar; colormap(redblue); caxis([-1*cmax cmax]);
set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
    'fontweight', 'normal','Position',[0.95 0.358 0.016 0.102])

set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
    'fontname','arial','fontsize',12,'Position',  [1.6 1.3*cmax 1])

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

hsub = subplot(4,2,8);
set(hsub, 'Position',[.57 .16 .34 .14])

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
    'fontweight', 'normal','Position',[0.95 0.168 0.016 0.102])

set(get(hcb,'title'),'string','(x10^{-7} m s^{-2})',...
    'fontname','arial','fontsize',12,'Position',  [1.6 1.3*cmax 1])



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
