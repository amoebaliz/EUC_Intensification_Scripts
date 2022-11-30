close all; clear all; clc

CO2 = ncread('fgco2_Omon_CESM1-BGC_esmHistorical_r1i1p1_185001-200512.nc','fgco2');
lat = ncread('fgco2_Omon_CESM1-BGC_esmHistorical_r1i1p1_185001-200512.nc','lat');
lon = ncread('fgco2_Omon_CESM1-BGC_esmHistorical_r1i1p1_185001-200512.nc','lon');

% indexing spatial bounds
 % I_lat = (find(lat(1,:)<=2.5 & lat(1,:)>=-2.5))';      
 I_lon = find(lon(:,1)>=209 & lon(:,1)<=269);
 I_lat = (find(lat(1,:)<=3.1 & lat(1,:)>=-3.1))';  
 
 slon = squeeze(lon(I_lon,1));
 slat = squeeze(lat(1,I_lat))';

 
 fgco2_trend = NaN(length(slon),length(slat));
 fgco2_plusminus = fgco2_trend;
 fgco2_sig = fgco2_trend;
 
 for j = 1:length(slon)
     for k = 1:length(slat)
     
     [fgco2_trend(j,k), fgco2_plusminus(j,k), fgco2_sig(j,k)] = trend_stat(squeeze(CO2(I_lon(j),I_lat(k),:)),95);
     
     end
 end
 
 fgco2_trend(fgco2_sig == 0) = NaN;
 
%%

close all
clc
subplot(2,1,1)

%%%%%%%%%%%%%% multiplied by 10^9
[a b] = contourf(slon,slat,10^9*squeeze(mean(CO2(I_lon,I_lat,:),3))',100);
set(b,'EdgeColor','none')

load redblue
colormap(redblue)
caxis([-1.5 1.5])
colorbar

% x axis
 xlim([210 270])
 xlabel ('Longitude', 'fontsize', 13, 'fontweight', 'bold') 
 xlabh = get(gca,'XLabel');
 set(xlabh,'Position',get(xlabh,'Position') - [0 1/3 0])
 set(gca,'Xtick',[150 180 210 240 270],'XTickLabel',{'150 E', '180', '150 W', '120 W', '90 W'},...
    'fontsize', 12, 'fontweight', 'bold')

% y axis
 ylim([-3 3])
 ylabel ('Latitude', 'fontsize', 13, 'fontweight', 'bold') 
 ylabh = get(gca,'XLabel');
  set(gca,'Ytick',[-2 0 2],'YTickLabel',{'-2 S', 'Eq', '1 S'},...
    'fontsize', 13, 'fontweight', 'bold')

subplot(2,1,2)

[a b] = contourf(slon,slat,10^12*fgco2_trend',100);
set(b,'EdgeColor','none')

hold on

h = contour(slon,slat,10^9*squeeze(mean(CO2(I_lon,I_lat,:),3))',10,'k');
% set(h,'EdgeColor','k')

colorbar
caxis([-.3 0.3])

 xlim([210 270])
 xlabel ('Longitude', 'fontsize', 13, 'fontweight', 'bold') 
 xlabh = get(gca,'XLabel');
 set(xlabh,'Position',get(xlabh,'Position') - [0 1/3 0])
 set(gca,'Xtick',[150 180 210 240 270],'XTickLabel',{'150 E', '180', '150 W', '120 W', '90 W'},...
    'fontsize', 12, 'fontweight', 'bold')

% y axis
 ylim([-3 3])
 ylabel ('Latitude', 'fontsize', 13, 'fontweight', 'bold') 
 ylabh = get(gca,'XLabel');
  set(gca,'Ytick',[-2 0 2],'YTickLabel',{'-2 S', 'Eq', '1 S'},...
    'fontsize', 13, 'fontweight', 'bold')

