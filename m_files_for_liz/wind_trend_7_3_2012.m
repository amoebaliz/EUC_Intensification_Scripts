% wind trends stuff

close all 
clear all
clc
addpath('/data1/mcode/')
addpath('/data1/mcode/seawater_ver3_3/')
load redblue.mat

ncid = netcdf.open ('SODA_2.2.6_Extended_Pac_taux.cdf','NC_NOWRITE');
ncid2 = netcdf.open ('SODA_2.2.6_Extended_Pac_tauy.cdf','NC_NOWRITE');

% get variable IDs for SEC
    varid_lat = netcdf.inqVarID(ncid,'LAT72_231');
    varid_lon = netcdf.inqVarID(ncid,'LON201_600');
    varid_taux = netcdf.inqVarID(ncid,'TAUX_ENS_MN');
    varid_tauy = netcdf.inqVarID(ncid2,'TAUY_ENS_MN');
    
% get variables (general)
    lat = netcdf.getVar(ncid,varid_lat);
    lon = netcdf.getVar(ncid,varid_lon);
    taux = netcdf.getVar(ncid,varid_taux);
    tauy = netcdf.getVar(ncid2,varid_tauy);

    taux(taux<-10000) = NaN;
    tauy(tauy<-10000) = NaN;

crl = windstresscurl(lon,lat,permute(taux,[3,2,1]),permute(tauy,[3,2,1]));
crl = permute(crl,[3,2,1]);

    
% calculate wind speed
wnd_speed = sqrt(taux.^2+tauy.^2);

avg_taux = mean(taux(:,:,:),3);
avg_tauy = mean(taux(:,:,:),3);

% calculate wind speed trend
wnd_trend = NaN(length(lon),length(lat));
wnd_plusmius = wnd_trend;
wnd_sig = wnd_trend;

crl_trend = wnd_trend;
crl_plusminus = wnd_trend;
crl_sig = wnd_trend;

taux_trend = wnd_trend;
taux_plusminus = taux_trend;
taux_sig = taux_trend;

for j = 1:length(lon)
    for k = 1:length(lat)
 
        [wnd_trend(j,k), wnd_plusminus(j,k), wnd_sig(j,k)] =trend_stat(squeeze(wnd_speed(j,k,:)),95);
        [crl_trend(j,k), crl_plusminus(j,k), crl_sig(j,k)] =trend_stat(squeeze(crl(j,k,:)),95);
        [taux_trend(j,k), taux_plusminus(j,k), taux_sig(j,k)] = trend_stat(squeeze(taux(j,k,:)),95);
    end
end



for j = 1:length(lon)
    for k = 1:length(lat)
 
        
    end
end
%%
close all

figure(1)
[~,h]=contourf(lon,lat, wnd_trend',100);
set(h,'EdgeColor','none')
colormap(redblue);
hold on

g = quiver(lon, lat, squeeze(mean(taux,3))',squeeze(mean(tauy,3))','k');

[xcoast,ycoast]=getcoast;

plot(xcoast,ycoast,'k')
%%
close all
figure(2)
[~,h]=contourf(lon,lat, (10^-8)*crl_trend',100);
set(h,'EdgeColor','none')
colormap(redblue);
hold on
[cs,q] = contour(lon,lat, taux_trend',-1*(1/10000)*[1 2 3 4],'-.b','LineWidth',2);
[cs2, qq] = contour(lon,lat, taux_trend',[0, 0, 0],':k', 'LineWidth',3);
[cs3, qqq] = contour(lon,lat, taux_trend',(1/10000)*[1 2 3 4],'k', 'LineWidth',2);
clabel(cs, q, 'FontSize', 10, 'Color', 'b');
clabel(cs3, qqq, 'FontSize', 10, 'Color', 'k');

[xcoast,ycoast]=getcoast;

caxis([-.000000001 0.000000001])

plot(xcoast,ycoast,'k')