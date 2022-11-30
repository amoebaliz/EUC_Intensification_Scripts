% SODA Scripts for SEC transports and taux 
%close all
clear all
clc

lat = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
lon = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
time = ncread('SODA_2.2.6_Trop_Pac_u.cdf','TIME');
depth = ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');
u = ncread('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN');

% eliminate missing value fill value (-9.99e+33)
u(u>1000) = NaN;

% Indexing spatial bounds
 I_lat = find(lat<=2 & lat>=-2);     
 I_lon = find(lon>=149 &lon<=271);
 I_dep = find(depth==depth);
 time0 = 61;
 
 slat = lat(I_lat);
 slon = lon(I_lon);
 sdep = depth(I_dep);

 nt = length(time);
 ndep = length(sdep);
 nlon = length(slon);
 nlat= length(slat);
 
%% zonal velocities averaged over latitudes (and time)
u_euc_eq = squeeze(nanmean(u(I_lon,I_lat,I_dep,time0:end),2));
u_euc_avg = squeeze(nanmean(nanmean(u(I_lon,I_lat,I_dep,time0:end),4),2));

%% statistical analyses
euc_trend = -500*ones(ndep,nlon);
euc_plusminus = euc_trend;
euc_sig = euc_trend;

%  i_mon_mam = zeros(35,3);
%  i_mon_jja = zeros(35,3);
% 
%  for i_nyr = 1:1656/12;
%                     
%     i_mon_mam(i_nyr,:) = (3:5)+12*(i_nyr-1);
%     i_mon_jja(i_nyr,:) = (6:8)+12*(i_nyr-1);
%     
%  end
% i_mon_mam = sort(i_mon_mam(:));
% i_mon_jja = sort(i_mon_jja(:));
% 
% t = i_mon_mam;

%close all
for j = 1: nlon
    for k = 1:ndep
        
        [euc_trend(k,j),euc_plusminus(k,j),euc_sig(k,j)]=trend_stat(squeeze(u_euc_eq(j,k,:)),99);
   
    end
end

euc_trend(euc_sig<1) = NaN;


save 'figure3_new.mat' 'euc_trend' 'u_euc_avg' 'slat' 'slon' 'sdep'

%% plotting data
close all
clc

figure('units','normal','outerposition',[.1 .5 .3 .6],'Color',[1 1 1])

load redblue
colormap(redblue)



[a b] = contourf(slon,sdep,12*100*euc_trend,500);
set(b,'EdgeColor','none')

cmin = -0.3;
cmax = 0.3;

caxis([cmin cmax])

hold on

[q,h] = contour(slon,sdep,u_euc_avg',[-.40,-.20],'--b','LineWidth',2);
[qq, hh] = contour(slon,sdep,u_euc_avg',[0, 0, 0],':k', 'LineWidth',2);
[qqq,hhh] = contour(slon,sdep,u_euc_avg',[.20,.40,.60,.80],'k', 'LineWidth',2);

set(gca,'YDir','reverse')
%%
ylim([0 500])
xlim([150 270])

set(gca,'TickDir','out')

ylabel ('Depth (m)', 'fontsize', 16, 'fontweight', 'bold') 
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [4 0 0])
set(gca,'YTick',[0 100 200 300 400 500])



xlabel ('Longitude', 'fontsize', 16, 'fontweight', 'bold') 
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [0 30 0])
%%
set(gca,'Xtick',[150 180 210 240 270],'XTickLabel',{'150 E', '180', '150 W', '120 W', '90 W'},...
    'fontsize', 13, 'fontweight', 'bold')

%clabel(q,h, 'manual','fontsize', 13,'color', 'blue', 'fontweight', 'bold','rotation',0)
%clabel(qq,hh, 'manual','fontsize', 13,'rotation',0)
%clabel(qqq,hhh, 'manual','fontsize', 13,'fontweight', 'bold','rotation',0)

%%
% set(gcf, 'Renderer', 'painters')
set(gca,'box','off')
isholdonque = ishold; 
hold on
ax = axis;
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',0.5)
plot(ax(1:2),1*[1,1],'k','linewidth',0.5)
if isholdonque == 0
hold off
end 

pos = get(gca,'OuterPosition');
set(gca,'OuterPosition',[pos(1) (pos(2)+.05) pos(3) (pos(4) - .1)])
%%
hcb = colorbar;
set(hcb,'YTick',[cmin, cmin/2, 0, cmax/2, cmax],'fontsize', 10,'fontweight', 'normal')
