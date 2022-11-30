% SODA Scripts for SEC transports and taux
 
 close all
 clear all
 
    % opening the 10 Lat to -10 lat info for current/ windstress information information
    ncid = netcdf.open ('u_eqpac_SODA_extended_lat10.nc','NC_NOWRITE');
    ncid_2 = netcdf.open ('taux_eqpac_SODA_extended_lat10.nc','NC_NOWRITE');    
    % for looking at EUC too
    ncid_3 = netcdf.open ('u_eqpac_SODA_extended.nc','NC_NOWRITE'); 
%%    
    % get variable IDs for SEC
    varid_lat_1 = netcdf.inqVarID(ncid,'lat');
    varid_lon = netcdf.inqVarID(ncid,'lon');
    varid_time = netcdf.inqVarID(ncid,'time');
    varid_u_sec = netcdf.inqVarID(ncid,'u');
    varid_depth_sec = netcdf.inqVarID(ncid,'depth');
    
    % get variable IDs for taux
    varid_taux = netcdf.inqVarID(ncid_2,'taux');
    
    % get variable IDs for EUC
    varid_u_euc = netcdf.inqVarID(ncid,'u');
    varid_depth_euc = netcdf.inqVarID(ncid,'depth');
    varid_lat_2 = netcdf.inqVarID(ncid,'lat');
    
    % get variables (general)
    time = netcdf.getVar(ncid,varid_time);
    sec_lat = netcdf.getVar(ncid,varid_lat_1);
    euc_lat = netcdf.getVar(ncid_3,varid_lat_2);
    lon = netcdf.getVar(ncid,varid_lon);
    
    % variables for SEC (lon x lat x depth x time)
    u_sec = netcdf.getVar(ncid,varid_u_sec);
    depth_sec = netcdf.getVar(ncid,varid_depth_sec);
    
    % variables for taux (lon x lat x time)
    taux = netcdf.getVar(ncid_2,varid_taux);
    
    % variables for EUC (lon x lat x depth x time)
    
    u_euc = netcdf.getVar(ncid_3,varid_u_euc);
    depth_euc = netcdf.getVar(ncid_3,varid_depth_euc);
 
 % eliminate missing value fill value (-9.99e+33)
 u_sec(u_sec<-10000) = NaN;
 taux(taux<-10000) = NaN;
 u_euc(u_euc<-10000) = NaN;
 
 %% Indexing spatial bounds
 
 depth_sec_2 = 15.0700;
 
 I_lat_1 = find(sec_lat<=7 & sec_lat>=-7);      % For SEC & taux
 I_lat_2 = find(euc_lat<=0.5 & euc_lat>=-0.5);      % For EUC
 I_lon = find(lon>=150 &lon<=270);
 I_dep = find(depth_euc>=50 & depth_euc<=300);
 
 slat_1 = sec_lat(I_lat_1);
 slat_2 = euc_lat(I_lat_2);
 slon = lon(I_lon);

 nt = length(time);
 nlon = length(slon);
 nlat_1 = length(slat_1);
 nlat_2 = length(slat_2);

%% Putting time in datenum form

time_2 = zeros(1,length(time));
 
 timen = time-.5;
 for j = 1:length(time)
    time_2(j) = addtodate(datenum(1960,1,1),double(timen(j)),'month');
 end
% fig1 = figure('Color',[1 1 1]);
% ss=get(0,'ScreenSize');
% set(fig1,'Position',[ss(1)+5 ss(2) ss(3)-250 ss(4)-120]) %%

u_euc2 = squeeze(nanmean(u_euc(:,I_lat_2,:,:),2));
%%
u_euc3 = -500*ones(size(u_euc2,1),size(u_euc2,2),size(u_euc2,3)/12);
%%

for j = 1: (size(u_euc2,3)/12)
    u_euc3(:,:,j) = nanmean(u_euc2(:,:,((j-1)*12+1):(j*12)),3);
end

u_euc3 = u_euc3(:,:,:)*100;
%u_euc4 = squeeze(nanmean(nanmean(u_euc(:,I_lat_2,:,:),2),4));


%%
 u_euc55 = squeeze(squeeze(nanmean(nanmean(u_euc(:,I_lat_2,:,:),2),4)));
% save('average_euc.mat', 'u_euc55', 'depth_euc', 'lon')


%%
% initializing variables for statistical analysis
euc_trend = -500*ones(length(depth_euc(3:18)),length(I_lon));euc_plusminus = euc_trend;euc_sig = euc_trend;

euc_nu = u_euc2(I_lon,:,:);

close all
for j = 1: length(I_lon)
    for k = 3:18
        
    [euc_trend(k-2,j),euc_plusminus(k-2,j),euc_sig(k-2,j)]=trend_stat((squeeze(squeeze(euc_nu(j,k,:)))),95);
    end
end

close all
euc_trend(euc_sig<1) = NaN;

save('filled_contours.mat', 'euc_trend')
%%
clc
close all
fifty_mean = squeeze(mean(u_euc3(:,:,:),3));
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

[a b] = contourf(slon,depth_euc(3:18),12*100*euc_trend,300);
set(b,'EdgeColor','none')
%t = pcolor(repmat(slon',[16,1]),repmat(depth_euc(3:18),[1,240]),double(100*euc_trend));

hold on
% fifty_mean = u_euc4(I_lon,:)*100; (same thing apparently)
q = contour(slon,depth_euc(3:18),fifty_mean(I_lon,3:18)',[-40,-20],'-.b','LineWidth',2);
qq = contour(slon,depth_euc(3:18),fifty_mean(I_lon,3:18)',[0, 0, 0],':k', 'LineWidth',2);
qqq = contour(slon,depth_euc(3:18),fifty_mean(I_lon,3:18)',[20,40,60,80],'k', 'LineWidth',2);
set(gca,'YDir','reverse')
%set(t, 'EdgeColor', 'none')
caxis([-.3 0.3])
xlim([150 270])
ylim([depth_euc(3) depth_euc(18)])

colorbar


%%
close all
I_time_1 = (75*12+1):(79*12);
I_time_2 = (108*12+1):(110*12);

sec_strong = squeeze(squeeze(nanmean(nanmean(u_euc(I_lon,I_lat_2,:,I_time_1),4),2)));
sec_weak = squeeze(squeeze(nanmean(nanmean(u_euc(I_lon,I_lat_2,:,I_time_2),4),2)));
% plot a strong SEC year
subplot(2,1,1)
[a b] = contourf(slon,depth_euc(1:10),sec_strong(:,1:10)',300);
%[a b] = contourf(slon,depth_euc,sec_strong(:,:)',300);
set(b,'EdgeColor','none')
set(gca,'YDir','reverse')
caxis([-.5 0.5])
colorbar
title('1946:1949')
%plot a weak SEC year
subplot(2,1,2)
[a b] = contourf(slon,depth_euc(1:10),sec_weak(:,1:10)',300);
%[a b] = contourf(slon,depth_euc,sec_weak(:,:)',300);
set(b, 'EdgeColor','none')
set(gca,'YDir','reverse')
caxis([-.5 0.5])
colorbar
title('1979:1981')

%%
close all



max_u = squeeze(squeeze(squeeze(max(max(max(u_euc(I_lon,I_lat_2,I_dep,:)))))));
% max_u_2 = squeeze(squeeze(squeeze(nanmean(max(max(u(I_lon,I_lat,I_dep,:),[],2),[],3)))));

 titles = {'monthly data','1-yr mean','7-yr mean','10-yr mean'};
 win_i = [0,1,7,10]; % number of years in running mean
 
 time_3 = time_2;
 
 close all
  for j = 1:4
 
      % get running mean  
      win = win_i(j)*12+1;
      avgmax_u = runmean(max_u,win);timeser = avgmax_u;
      
      % calculate trend
      trend= timeser((1+floor(win/2)):(end-floor(win/2)))-detrend(timeser((1+floor(win/2)):(end-floor(win/2))));
      
     subplot(2,2,j)
     
     % plot data or running mean
     plot(time_3,avgmax_u)
     hold on
     % plot trendline
     plot(time_3((1+floor(win/2)):(end-floor(win/2))),trend,'-.r','LineWidth',2)
     
     %plot details
     legend(char(titles(j)),'Location', 'northwest')
     title((trend(2)-trend(1))*12*100); % this calculates the trend per century if the data are daily
     datetick('x','yyyy','keeplimits')
    
 end