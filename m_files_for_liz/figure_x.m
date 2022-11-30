close all; clc

load redblue
load ududx_trend.mat
load vdudy_trend.mat
load wdudz_trend.mat
load zpgf_trend.mat


trend_sum = ududx_trend(2:end-1,:) + vdudy_trend(2:end-1,2:end-1) + ...
    wdudz_trend(:,2:end-1) + zpgf_trend(2:end-1,:);

%% ududx trend plot
subplot(7,1,1)
colormap(redblue)

ududx_trend(ududx_sig == 0) = NaN;
[c, h] = contourf(slon(2:end-1),depth2',12*100*ududx_trend,500);
set(gca, 'YDir', 'reverse')
set(h,'EdgeColor','none')

hold on

[q, h] = contour(slon(2:end-1),depth2',avg_ududx',[0 0 0],'k','LineWidth',2);
[qq, hh] = contour(slon(2:end-1),depth2',avg_ududx',(.5:.5:2)*10^(-7),'k');
[qqq, hhh] = contour(slon(2:end-1),depth2',avg_ududx',(.5:.5:2)*-10^(-7),':k');

cmax = 5*10^(-4); cmin = -1*cmax; caxis([cmin cmax])

title('ududx')
colorbar

%% vdudy trend plot
subplot(7,1,2)
colormap(redblue)

vdudy_trend(vdudy_sig == 0) = NaN;
[c, h] = contourf(slon,depth2',12*100*vdudy_trend,500);
set(gca, 'YDir', 'reverse')
set(h,'EdgeColor','none')

hold on

[q, h] = contour(slon,depth2',avg_vdudy',[0 0 0],'k','LineWidth',2);
[qq, hh] = contour(slon,depth2',avg_vdudy',(.5:.5:2)*10^(-7),'k');
[qqq, hhh] = contour(slon,depth2',avg_vdudy',(.5:.5:2)*-10^(-7),':k');

cmax = 5*10^(-4); cmin = -1*cmax; caxis([cmin cmax])

title('vdudy')
colorbar

%% wdudz trend plot
subplot(7,1,3)
colormap(redblue)

wdudz_trend(wdudz_sig == 0) = NaN;
[c, h] = contourf(slon,depth2(2:end-1)',12*100*wdudz_trend,500);
set(gca, 'YDir', 'reverse')
set(h,'EdgeColor','none')

hold on

[q, h] = contour(slon,depth2(2:end-1)',avg_wdudz',[0 0 0],'k','LineWidth',2);
[qq, hh] = contour(slon,depth2(2:end-1)',avg_wdudz',(.5:.5:2)*10^(-7),'k');
[qqq, hhh] = contour(slon,depth2(2:end-1)',avg_wdudz',(.5:.5:2)*-10^(-7),':k');

cmax = 5*10^(-4); cmin = -1*cmax; caxis([cmin cmax])

title('wdudz')
colorbar


%% zpgf trend plot
subplot(7,1,4)
colormap(redblue)

zpgf_trend(zpgf_sig == 0) = NaN;
[c, h] = contourf(slon(2:end-1),depth2',12*100*zpgf_trend,500);
set(gca, 'YDir', 'reverse')
set(h,'EdgeColor','none')

hold on

[q, h] = contour(slon(2:end-1),depth2',avg_zpgf',[0 0 0],'k','LineWidth',2);
[qq, hh] = contour(slon(2:end-1),depth2',avg_zpgf',(.5:.5:2)*10^(-7),'k');
[qqq, hhh] = contour(slon(2:end-1),depth2',avg_zpgf',(.5:.5:2)*-10^(-7),':k');

cmax = 5*10^(-4); cmin = -1*cmax; caxis([cmin cmax])

title('zpgf')
colorbar

%% frich trend plot
subplot(7,1,5)

%% fricv trend plot
subplot(7,1,6)

%% sum of trends
subplot(7,1,7)
colormap(redblue)

[c, h] = contourf(slon(2:end-1),depth2(2:end-1)',12*100*trend_sum,500);
set(gca, 'YDir', 'reverse')
set(h,'EdgeColor','none')

cmax = 5*10^(-4); cmin = -1*cmax; caxis([cmin cmax])

title('sum of trends')
colorbar