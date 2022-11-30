close all
clear all
clc

lat = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
lon = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
time = ncread('SODA_2.2.6_Trop_Pac_u.cdf','TIME');
depth = ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');
u = ncread('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN');

 % eliminate missing value fill value (-9.99e+33)
u(u>1000) = NaN;


 % indexing spatial bounds
 I_lat = find(lat<=2 & lat>=-2);      % For EUC
 
 I_dep = find(depth>=10 & depth<=300);
 
 %actual values for EUC
nlon = length(lon); 

trend_sens=NaN(nlon,nlon);
u_sig = NaN(nlon,nlon);

for j = 1:nlon
    for k = j+1:nlon
        
    I_lon = find(lon>=lon(j) & lon<=lon(k));
        
%% find max zonal velocity 
    max_u = squeeze(max(max(max(u(I_lon, I_lat,I_dep,:)))));

      
%% Calculate the trend

    [u_trend,u_plusminus,u_sig(j,k)]=trend_stat(max_u,99);
    
    trend_sens(j,k)=u_trend*12*100;
    
    end
end
%%
close all
[~,h]=contourf(lon,lon,trend_sens',100);
set(h,'EdgeColor','none')
colorbar
caxis([0 .25])
%set(gca,'YDir','reverse')
set(gca,'FontSize',16)
xlabel('West Longitudinal Limit'); ylabel('East Longitudinal Limit'); 
title('Sensitivity Analysis of Impact of Sampled Longitudinal Range on EUC Trend')