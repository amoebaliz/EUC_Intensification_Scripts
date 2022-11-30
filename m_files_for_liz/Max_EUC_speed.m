close all
clear all
clc

lat = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
lon = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
time = ncread('SODA_2.2.6_Trop_Pac_u.cdf','TIME');
depth = ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');
u = ncread('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN');

 % eliminate missing value fill value (-1e+34)
u(u<-10000) = NaN;

 % indexing spatial bounds
 I_lat = find(lat<=2 & lat>=-2);      % For EUC
 I_lon = find(lon>=210 &lon<=270);
 I_dep = find(depth>=10 & depth<=300);
 time0 = 61;
 
 %actual values for EUC
 slat = lat(I_lat);
 slon = lon(I_lon);

 nt = length(time(time0:end));
 nlon = length(slon);
 nlat = length(slat);
 
 %% Putting time in datenum form

time_2 = zeros(1,length(time(time0:end)));
 
 for j = time0:length(time)
    time_2(j-time0+1) = addtodate(datenum(1865,12,15),double(time(j)),'month');
 end

%% find max zonal velocity 
max_u = squeeze(max(max(max(u(I_lon, I_lat,I_dep,time0:end)))));

F = filtrage2(max_u,'low',9,7*12);
F2 = filtrage2(max_u,'low',9,10*12);

%save('max_EUC_runmean.mat', 'avgmax_u')
      
%% Calculate the trend

[u_trend,u_plusminus,u_sig]=trend_stat(max_u,99);

time_3 = datevec(time_2);

save ('max.mat', 'time_2', 'u_trend', 'F' , 'max_u')

figure
subplot(2,1,1)
h = plot(time_2,max_u);
hold on
z = plot(time_2,F,'k','LineWidth',2);
datetick('x', 10,'keepticks')
xlabel('Year')
ylabel('Meters Per Second (m/s)')
title('PER = 7*12')

subplot(2,1,2)
h = plot(time_2,max_u);
hold on
z = plot(time_2,F2,'k','LineWidth',2);
datetick('x', 10,'keepticks')
xlabel('Year')
ylabel('Meters Per Second (m/s)')
title('PER = 10*12')