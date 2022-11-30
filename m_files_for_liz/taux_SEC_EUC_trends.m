% SODA Scripts for transports and taux
 
 close all
 %clear all
 
lat = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
lon = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
time = ncread('SODA_2.2.6_Trop_Pac_u.cdf','TIME');
depth = ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');
%%
u = ncread('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN');
taux = ncread('SODA_2.2.6_Trop_Pac_taux.cdf','TAUX_ENS_MN');

% eliminate missing value fill value (-9.99e+33)
u(u>1000) = NaN;
taux(taux>10000) = NaN;

 % Indexing spatial bounds

 I_lat = find(lat<=0.5 & lat>=-0.5);     
 I_lon = find(lon>=150 &lon<=270);
 I_dep = find(depth>=10 & depth<=300);
 %I_dep = find(depth==depth);
 time0 = 61;
 
 slat = lat(I_lat);
 slon = lon(I_lon);

 nt = length(time);
 nlon = length(slon);
 nlat= length(slat);
 
%% Building Transport Areal Grids
 
% calculating area of grid box; depth dimension:

depth_height = zeros(length(depth),1);

for d = 1:length(depth)
    
    depth_height(d) = 2*(depth(d)-sum(depth_height));
    
end

% latitude dimension:
% 
% lat1 = slat+.25;
% lat2 = slat-.25;

%Convert to radians

% latrad1 = lat1*pi/180;
% latrad2 = lat2*pi/180;
% 
% raddis = acos(sin(latrad2).*sin(latrad1) + cos(latrad2).*cos(latrad1)*1);
% nautdis = raddis * 3437.74677;
% stdiskm = nautdis * 1.852;
% stdism = stdiskm * 1000;

lat_dim = 0.5*110.575*1000;       %in meters

% calculate grid of areas (square meters)

grid_area = repmat(depth_height*lat_dim,[1,length(I_lat),length(I_lon),length(time)]);
grid_area = permute(grid_area,[3 2 1 4]);

%% Calculate transports and average wind stress

%current_transport = grid_area.*u(I_lon,I_lat,:,:);

% IN SVERDRUPS
% SEC
%SEC_transport = squeeze(nansum(current_transport(:,:,1,:),2))/(10^6);

SEC_velocity = squeeze(mean(u(I_lon,I_lat,1,time0:end),2));

% EUC
%EUC_transport = squeeze(nansum(nansum(current_transport(:,:,2:end,:),3),2))/(10^6);
EUC_velocity = squeeze(max(squeeze(mean(u(I_lon,I_lat,I_dep,time0:end),2)),[],2));

% TAUX
avg_taux = squeeze(mean(taux(I_lon,I_lat,time0:end),2));

%% Calculating Statistical Trends

% initializing variables for statistical analysis
sec_trend = NaN(length(I_lon),1); sec_plusminus = sec_trend; sec_sig = sec_trend;
euc_trend = sec_trend; euc_plusminus = euc_trend;euc_sig = euc_trend;
taux_trend = sec_trend; taux_plusminus = taux_trend; taux_sig = taux_trend;

for j = 1: length(I_lon)
    [sec_trend(j),sec_plusminus(j),sec_sig(j)]=trend_stat(SEC_velocity(j,:)',99);
    [euc_trend(j),euc_plusminus(j),euc_sig(j)]=trend_stat(EUC_velocity(j,:)',99);
    [taux_trend(j),taux_plusminus(j),taux_sig(j)]=trend_stat(avg_taux(j,:)',99);
end

save('trend_data.mat', 'euc_trend', 'sec_trend', 'taux_trend', 'sec_plusminus',...
    'euc_plusminus', 'taux_plusminus','slon')

%% Create Plots
% close all
% per century conversion
k = 12*100;
figure
subplot(3,1,1)
errorbar(slon,taux_trend*k,taux_plusminus*k,'o', 'MarkerFace','k')
v=axis;
hold on
plot(v(1:2),[0 0],'--r','LineWidth',2)
subplot(3,1,2)
errorbar(slon,sec_trend*k,sec_plusminus*k,'o', 'MarkerFace','k')
v=axis;
hold on
plot(v(1:2),[0 0],'--r','LineWidth',2)
subplot(3,1,3)
errorbar(slon,euc_trend*k,euc_plusminus*k,'o', 'MarkerFace','k')
v=axis;
hold on
plot(v(1:2),[0 0],'--r','LineWidth',2)



 
 