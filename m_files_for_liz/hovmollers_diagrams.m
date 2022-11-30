close all
clear all
clc

time0 = 61;

lon = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
lat = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
depth = ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');
time = ncread('SODA_2.2.6_Trop_Pac_u.cdf','TIME');
time = time(time0:end);

taux = ncread('SODA_2.2.6_Trop_Pac_taux.cdf','TAUX_ENS_MN');
tuax = taux(:,:,time0:end);

ssh = ncread('SODA_2.2.6_Trop_Pac_ssh.cdf','SSH_ENS_MN');
ssh = ssh(:,:,time0:end);

u = ncread('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN');
u = u(:,:,:,time0:end);

temp =  ncread('SODA_2.2.6_Trop_Pac_temp.cdf','TEMP_ENS_MN');
temp = temp(:,:,:,time0:end);

% eliminate missing value fill value (-9.99e+33)
taux(taux<-10000) = NaN;
ssh(ssh<-10000) = NaN;
u(u<-10000) = NaN;
temp(temp<-10000) = NaN;
%
% Indexing spatial bounds
 I_lat = find(lat<=2 & lat>=-2);     
 I_lon = find(lon>=149 &lon<=271);
 I_dep = find(depth >10 & depth<300);

 slat = lat(I_lat);
 slon = lon(I_lon);
 sdep = depth(I_dep);

 nt = length(time);
 ndep = length(sdep);
 nlon = length(slon);
 nlat= length(slat);
 
%%
% calculating area of grid box; depth dimension:
depth_height = zeros(length(depth),1);

for d = 1:length(depth)
    depth_height(d) = 2*(depth(d)-sum(depth_height));
end

lat_dim = 0.5*110.575*1000;       %in meters

% calculate grid of areas (square meters)

grid_area = repmat(depth_height*lat_dim,[1,nlat,nlon,nt]);
grid_area = permute(grid_area,[3 2 1 4]);

u_transport = u(I_lon,I_lat,:,:).*grid_area;

max_euc = squeeze(max(max(u(I_lon,I_lat,I_dep,:),[],2),[],3));

%%
% climatologies

 [taux_clim, ~] = climanom(permute(mean(taux(I_lon,I_lat,:),2),[3,2,1]));
 [temp_clim, ~] = climanom(permute(mean(squeeze(temp(I_lon,I_lat,1,:)),2),[3,2,1]));
 [ssh_clim, ~] = climanom(permute(mean(ssh(I_lon,I_lat,:),2),[3,2,1]));
 [sec_clim, ~] = climanom(permute(mean(squeeze(u(I_lon,I_lat,1,:)),2),[3,2,1]));
 
 [euc_clim, ~] = climanom(permute(mean(squeeze(sum(u_transport(:,:,:,:),3)),2),[3,2,1]));
 [max_euc_clim,~] = climanom(permute(max_euc,[2 3 1]));
 

%%
% trends
I_month = 1:12:(nt-11);

taux_trend = zeros(12,nlon); taux_sig = zeros(12,nlon);
temp_trend = zeros(12,nlon); temp_sig = zeros(12, nlon);
ssh_trend = zeros(12,nlon); ssh_sig = zeros(12,nlon);
sec_trend = zeros(12,nlon); sec_sig = zeros(12,nlon);
euc_trend = zeros(12,nlon); euc_sig = zeros(12,nlon);
max_euc_trend = zeros(12,nlon); max_euc_sig = zeros(12,nlon);

for j = 1: 12
    for k = 1:nlon
        
        [taux_trend(j,k),~,taux_sig(j,k)]=trend_stat(squeeze(mean(taux(k,I_lat,(I_month + (j-1))),2)),95);
        [temp_trend(j,k),~,temp_sig(j,k)]=trend_stat(squeeze(mean(temp(k,I_lat,1,(I_month + (j-1))),2)),95);
        [ssh_trend(j,k),~,ssh_sig(j,k)]=trend_stat(squeeze(mean(ssh(k,I_lat,(I_month + (j-1))),2)),95);       
        [sec_trend(j,k),~,sec_sig(j,k)]=trend_stat(squeeze(mean(u(k,I_lat,1,(I_month + (j-1))),2)),95);
        [euc_trend(j,k),~,euc_sig(j,k)]=trend_stat(squeeze(sum(mean(u_transport(k,:,:,(I_month + (j-1))),2),3)),95);
        [max_euc_trend(j,k),~,max_euc_sig(j,k)] = trend_stat(squeeze(max_euc(k,(I_month + (j-1))))',95);
   
    end
end
%%


%figure

 taux_trend(taux_sig<1) = NaN;
 temp_trend(temp_sig<1) = NaN;
 ssh_trend(ssh_sig<1) = NaN;
 sec_trend(sec_sig<1) = NaN;
 euc_trend(euc_sig<1) = NaN;
 max_euc_trend(max_euc_sig<1) = NaN;

%% save matrix of data
save hovs95.mat taux_clim temp_clim ssh_clim sec_clim euc_clim  max_euc_clim... 
    taux_trend temp_trend ssh_trend sec_trend euc_trend euc_trend max_euc_trend slon