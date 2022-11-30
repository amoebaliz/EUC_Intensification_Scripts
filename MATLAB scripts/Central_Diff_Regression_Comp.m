
% Checking the trend result when looking at central differencing of zonal
% velocity over time

close all 
clear all
clc

ncid = netcdf.open ('SODA_2.2.6_Trop_Pac_u.cdf','NC_NOWRITE');

% get variable IDs for SEC
    varid_lat = netcdf.inqVarID(ncid,'LAT142_161');
    varid_lon = netcdf.inqVarID(ncid,'LON241_560');
    varid_time = netcdf.inqVarID(ncid,'TIME');
    varid_u = netcdf.inqVarID(ncid,'U_ENS_MN');
    varid_depth = netcdf.inqVarID(ncid,'DEPTH1_20');

    % get variables (general)
    time = netcdf.getVar(ncid,varid_time);
    lat = netcdf.getVar(ncid,varid_lat);
    lon = netcdf.getVar(ncid,varid_lon);
    u = netcdf.getVar(ncid,varid_u);
    depth = netcdf.getVar(ncid,varid_depth);
    
 % eliminate missing value fill value (-9.99e+33)
 u(u<-10000) = NaN;
 
% Indexing spatial bounds

 I_lat = find(lat<=0.5 & lat>=-0.5);  
 %I_lon = find(lon>=150 &lon<=270);
 I_lon = 1:length(lon);
 
 nt = length(time);
 slon = lon(I_lon);
 
 u_eq = squeeze(mean(u(I_lon,I_lat,:,:),2));
 avg_u = 100*squeeze(nanmean(nanmean(u(I_lon,I_lat,:,:),2),4));
 
 %% calculating trend via central differencing - 1st method
% VERY VERY SLOW and gives same result as 2nd method
%  trend_val = NaN(length(I_lon),length(depth));
%  
%  for j = 1:length(I_lon)
%     for k = 1:length(depth)
%         cd = NaN(nt-2,1);
%         for t = 2:(length(time)-1)
%             cd(t-1) = (u_eq(j,k,t+1)-u_eq(j,k,t-1))/2; 
%         end
%         
%         trend_val(j,k) = sum(cd)/(nt-2);
%     end
% end
 
 %% calculating trend via central differencing - 2nd method
% central difference 
center_dif = (u_eq(:,:,3:end)-u_eq(:,:,1:end-2))/2;
% average slope
avg_dif = squeeze(sum(center_dif,3)/(nt-2));

 %% calculating trend via forward differencing - 2nd method
 
 % forward difference 
for_dif = (u_eq(:,:,2:end)-u_eq(:,:,1:(end-1)))/2;
 % average slope
avg_dif_2 =  squeeze(sum(for_dif,3)/(nt-2));


%% calculating trend via regression

% initializing variables for statistical analysis
euc_trend = NaN(length(depth),length(I_lon));
euc_plusminus = euc_trend; 
euc_sig = euc_trend;

for j = 1: length(I_lon)
    for k = 1:length(depth)
        
    [euc_trend(k,j),euc_plusminus(k,j),euc_sig(k,j)]=trend_stat((squeeze(squeeze(u_eq(j,k,:)))),99);
    end
end
 
%euc_trend(euc_sig<1) = NaN;
 %% Figure Generating
 
 %COLORMAP DEFINITION
colorStart = [0 .1 .9];           % Blue (negative end of spectrum)
colorCenter = [1 1 1];    % Grey (to be centered at 'caxis' 0)
colorEnd = [.9 0 0];             % Red  (positive end of spectrum)
num = 30;

cmap1 = zeros(num,3);
cmap2 = cmap1;

for j = 1:3
    
    cmap1(1:num,j)= linspace(colorStart(j), colorCenter(j),num);
    cmap2(1:num,j) = linspace(colorCenter(j), colorEnd(j),num);
end

fig1 = figure('Color',[1 1 1]);

set(fig1,'Position',[400 600 1100 600]) 
cmap = [cmap1(1:end-1,:);cmap2(:,:)];
colormap(cmap)
% 
% subplot(3,1,1)
% [a b] = contourf(slon,depth,trend_val',300);
% set(b,'EdgeColor','none')
% 
% hold on
% 
% q = contour(slon,depth,avg_u',[-40,-20],'-.b','LineWidth',2);
% qq = contour(slon,depth,avg_u',[0, 0, 0],':k', 'LineWidth',2);
% qqq = contour(slon,depth,avg_u',[20,40,60,80],'k', 'LineWidth',2);
% 
% 
% set(gca,'YDir','reverse')
% caxis([-.0003 0.0003])

subplot(3,1,1)
[a b] = contourf(slon,depth,12*100*avg_dif',300);
set(b,'EdgeColor','none')

hold on

q = contour(slon,depth,avg_u',[-40,-20],'-.b','LineWidth',2);
qq = contour(slon,depth,avg_u',[0, 0, 0],':k', 'LineWidth',2);
qqq = contour(slon,depth,avg_u',[20,40,60,80],'k', 'LineWidth',2);


set(gca,'YDir','reverse')
caxis([-.3 0.3])
title('Trend by Central Differencing u (m sec^{-1} century^{-1})')
colorbar

subplot(3,1,2)
[a b] = contourf(slon,depth,12*100*avg_dif_2',300);
set(b,'EdgeColor','none')

hold on

q = contour(slon,depth,avg_u',[-40,-20],'-.b','LineWidth',2);
qq = contour(slon,depth,avg_u',[0, 0, 0],':k', 'LineWidth',2);
qqq = contour(slon,depth,avg_u',[20,40,60,80],'k', 'LineWidth',2);


set(gca,'YDir','reverse')
caxis([-.3 0.3])
title('Trend by Forward Differencing u (m sec^{-1} century^{-1})')

subplot(3,1,3)

 
[a b] = contourf(slon,depth,12*100*euc_trend,300);
set(b,'EdgeColor','none')
hold on
q = contour(slon,depth,avg_u',[-40,-20],'-.b','LineWidth',2);
qq = contour(slon,depth,avg_u',[0, 0, 0],':k', 'LineWidth',2);
qqq = contour(slon,depth,avg_u',[20,40,60,80],'k', 'LineWidth',2); 
set(gca,'YDir','reverse')
caxis([-.3 0.3])
title('Trend by Regressing u (m sec^{-1} century^{-1})')
colorbar