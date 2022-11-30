% SODA Scripts for SEC transports and taux
 
 close all
 clear all
 
    ncid = netcdf.open ('u_eqpac_SODA_extended_lat10.nc','NC_NOWRITE');      % SEC
    ncid_2 = netcdf.open ('taux_eqpac_SODA_extended_lat10.nc','NC_NOWRITE'); % TAUX   
    ncid_3 = netcdf.open ('u_eqpac_SODA_extended.nc','NC_NOWRITE');          % EUC
    ncid_4 = netcdf.open ('ssh_eqpac_SODA_extended.nc','NC_NOWRITE');        % SSH
    
    % get variable IDs for SEC
    varid_lat_1 = netcdf.inqVarID(ncid,'lat');
    varid_lon = netcdf.inqVarID(ncid,'lon');
    varid_time = netcdf.inqVarID(ncid,'time');
    varid_u_sec = netcdf.inqVarID(ncid,'u');
    varid_depth_sec = netcdf.inqVarID(ncid,'depth');
   
    
    % get variable IDs for taux
    varid_taux = netcdf.inqVarID(ncid_2,'taux');
    
    % get variable IDs for EUC
    varid_u_euc = netcdf.inqVarID(ncid_3,'u');
    varid_depth_euc = netcdf.inqVarID(ncid_3,'depth');
    varid_lat_2 = netcdf.inqVarID(ncid_3,'lat');
    
    % get variable IDs for SSH
     varid_ssh = netcdf.inqVarID(ncid_4,'ssh');
    
    % get variables (general)
    time = netcdf.getVar(ncid,varid_time);
    sec_lat = netcdf.getVar(ncid,varid_lat_1);
    euc_lat = netcdf.getVar(ncid_3,varid_lat_2);
    lon = netcdf.getVar(ncid,varid_lon);
    
    % variables for SEC (lon x lat x depth x time)
    u_sec = netcdf.getVar(ncid,varid_u_sec);
    u_sec = squeeze(u_sec);
    depth_sec = netcdf.getVar(ncid,varid_depth_sec);
    
    % variables for taux (lon x lat x time)
    taux = netcdf.getVar(ncid_2,varid_taux);
    
    % variables for EUC (lon x lat x depth x time)
    
    u_euc = netcdf.getVar(ncid_3,varid_u_euc);
    depth_euc = netcdf.getVar(ncid_3,varid_depth_euc);
    
    % variables for SSH (lon x lat x depth x time)
    ssh = netcdf.getVar(ncid_4,varid_ssh);
 % eliminate missing value fill value (-9.99e+33)
 u_sec(u_sec<-10000) = NaN;
 taux(taux<-10000) = NaN;
 u_euc(u_euc<-10000) = NaN;
 

 %% Indexing spatial bounds
 
 % depth_sec_2 = 15.0700;
 
 I_lat_1 = find(sec_lat<=1 & sec_lat>=-1);      % For SEC & taux
 I_lat_2 = find(euc_lat<=2 & euc_lat>=-2);      % For EUC
 I_lat_3 = find(euc_lat<=1 & euc_lat>=-1);
 I_lon = find(lon>=160 &lon<=260);
 I_depth = find(depth_euc>=0 & depth_euc<=500);
 
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
 
 %% Building Transport Areal Grids
 
% calculating area of grid box:
% depth dimension
% height_1 = (depth_sec(1)+depth_sec_2)/2;    % SEC

height_2 = zeros(length(depth_euc),1);  % EUC
depth_int = diff([0;depth_euc]);

for j = 2:length(depth_euc)
    height_2(j-1) = 0.5*(depth_int(j) + depth_int(j-1));  
end
height_2(1) = height_2(1)+ 0.5*depth_int(1);
height_2(end) = depth_int(end);

 
% latitudinal dimension: always about 0.5 degrees apart, assume 110.574 km
% per degree lat
% ASSUMPTION: distant between latitudes is constant - not actually true

y_dim = 0.5*110.574*100; %in meters

% calculate grid of areas (square meters)

% SEC
% sec_grid_area_1 = repmat(height_1*y_dim,1,numel(find(I_lat_1)));
% sec_grid_area_2 = repmat(sec_grid_area_1,[1 1 length(lon)]);
% sec_grid_area_3 = permute(sec_grid_area_2,[3 2 1]);
% sec_grid_area_4 = repmat(sec_grid_area_3,[1 1 1 length(time)]);

% EUC
euc_grid_area_1 = repmat(height_2*y_dim,1,numel(find(I_lat_2)));
euc_grid_area_2 = repmat(euc_grid_area_1,[1 1 length(lon)]);
euc_grid_area_3 = permute(euc_grid_area_2,[3 2 1]);
euc_grid_area_4 = repmat(euc_grid_area_3,[1 1 1 length(time)]);


%% Calculate transports

% SEC
% sec_transport_1 = sec_grid_area_4.*u_sec(:,I_lat_1,:,:);
% looking at westward flow
% sec_transport_2 = sec_transport_1;
% sec_transport_2(sec_transport_1>0) = 0;
% sec_transport_2 = -1*sec_transport_2;
% sum along lat, depth, divide by 10^6
% sec_transport_3 = squeeze(squeeze(nansum(nansum(sec_transport_2,2),3)))/(10^6); 


% EUC transport
euc_transport_1 = euc_grid_area_4.*u_euc(:,I_lat_2,:,:);
% looking at eastward flow
euc_transport_2 = euc_transport_1;
euc_transport_2(euc_transport_1<0) = 0;
% sum along lat, depth, divide by 10^6
euc_transport_3 = squeeze(squeeze(nansum(nansum(euc_transport_2,2),3)))/(10^6); 

% SEC 2 (above 200m)
% sec2_transport=euc_transport_1(:,:,(depth_euc>200),:);
% looking at westward flow
% sec2_transport(sec2_transport>0)=0;
% look at absolute value of westward flow
% sec2_transport = -1*sec2_transport;


%% Calculate Average Windstress

avg_taux = squeeze(nanmean(taux(:,I_lat_1,:),2));
% avg_taux = -1*avg_taux;

%% Clculate Average Surface Curent Velocity Between +/-1 degrees north

avg_SEC = squeeze(nanmean(u_sec(:,I_lat_1,:),2));

%% Clculate Average Sea Surface Height Between +/-1 degrees north

avg_SSH = squeeze(nanmean(ssh(:,I_lat_3,:),2));

%% Calculate Average EUC Velocity Between +/-2 and 0-500 m depth

avg_EUC = squeeze(nanmean(nanmean(u_euc(:,I_lat_2,:,:),2),3));
%% Calculate max EUC velocity Between +/-2 and 0-500 m depth

max_EUC = squeeze(squeeze(max(max(u_euc(:,I_lat_2,:,:),[],2),[],3)));
%%

% Initialize Hovmoller matrizes: all will have 12 rows (months) and 320 
% columns (longitudes)
TAUX_hov = zeros(12,length(lon)); TAUX_sig = TAUX_hov;
SEC_hov = zeros(12,length(lon)); SEC_sig = SEC_hov;
u_EUC_hov = zeros(12,length(lon)); u_EUC_sig = u_EUC_hov;
SSH_hov = zeros(12,length(lon)); SSH_sig = SSH_hov;
transport_EUC_hov = zeros(12,length(lon)); transport_EUC_sig = transport_EUC_hov;

%% Calculating Statistical Trends ON MONTHS!

nyear = size(time,1)/12;
nmonths = 1:nyear;

for h = 1:12 % for each month
    imonths = (nmonths-1)*12+h;
    for j = 1:length(lon)
    % determine trend in a given month of time @ 99% confidence interval
    [TAUX_hov(h,j), plus_minus, TAUX_sig(h,j)] = trend_stat(avg_taux(j,imonths)',95);
    [SEC_hov(h,j), plus_minus, SEC_sig(h,j)] = trend_stat(avg_SEC(j,imonths)',95);
    [u_EUC_hov(h,j), plus_minus, u_EUC_sig(h,j)] = trend_stat(max_EUC(j,imonths)',95);
    [SSH_hov(h,j), plus_minus, SSH_sig(h,j)] = trend_stat(avg_SSH(j,imonths)',95);
    [transport_EUC_hov(h,j), plus_minus, transport_EUC_sig(h,j)] = trend_stat(euc_transport_3(j,imonths)',95);
    end
    
end
%% eliminate non-significant regions 
TAUX_hov(TAUX_sig==0)= nan;
SEC_hov(SEC_sig==0) = nan;
u_EUC_hov(u_EUC_sig==0) = nan;
SSH_hov(SSH_sig == 0) = nan;
transport_EUC_hov(transport_EUC_sig == 0) = nan;

%% SAVE in proper form!
k = 12*100;
HOV_data = zeros(3840,7);

HOV_data(:,1) = repmat(lon,12,1);

lon_data = repmat((1:12),320,1);
HOV_data(:,2) = lon_data(:);

TAUX_hov_2 = TAUX_hov';
HOV_data(:,3) = TAUX_hov_2(:);

SEC_hov_2 = SEC_hov';
HOV_data(:,4) = SEC_hov_2(:);

u_EUC_hov_2 = u_EUC_hov';
HOV_data(:,5) = u_EUC_hov_2(:);

SSH_hov_2 = SSH_hov';
HOV_data(:,6) = SSH_hov_2(:);

transport_EUC_hov_2 = transport_EUC_hov';
HOV_data(:,7) = transport_EUC_hov_2(:);

HOV_data(:,8:12)=HOV_data(:,3:7)*k;


save('HOV_data.mat', 'HOV_data')

%%

% testing hovmollers
close all

k = 12*100;
contourf(SSH_hov*k);
colorbar

%%

% % initializing variables for statistical analysis
% sec_trend = -500*ones(length(lon),1); sec_plusminus = sec_trend; sec_sig = sec_trend;
% euc_trend = sec_trend; euc_plusminus = euc_trend;euc_sig = euc_trend;
% %sec2_trend = sec_trend; sec2_plusminus = sec2_trend; sec2_sig = sec2_trend;
% taux_trend = sec_trend; taux_plusminus = taux_trend; taux_sig = taux_trend;
% 
% for j = 1: length(lon)
%     % [sec_trend(j),sec_plusminus(j),sec_sig(j)]=trend_stat(sec_transport_3(j,:)',99);
%     % [sec2_trend(j),sec2_plusminus(j),sec2_sig(j)]=trend_stat(sec2_transport(j,:)',99);
%     [euc_trend(j),euc_plusminus(j),euc_sig(j)]=trend_stat(euc_transport_3(j,:)',99);
%     [taux_trend(j),taux_plusminus(j),taux_sig(j)]=trend_stat(avg_taux(j,:)',99);
% end
% 
% close all
%

%save('SEC_trend_data.mat', 'sec_trend' , 'sec_plusminus' , 'sec_sig')
%save('EUC_trend_data.mat', 'euc_trend' , 'euc_plusminus' , 'euc_sig')
%save('TAUX_trend_data.mat', 'taux_trend' , 'taux_plusminus' , 'taux_sig')
%save('SODA_lon.mat','lon')


%% Create Plots
% close all
% % per century conversion
 k = 12*100;
% 
% subplot(3,1,1)
% errorbar(lon,taux_trend*k,taux_plusminus*k,'o', 'MarkerFace','k')
% v=axis;
% hold on
% plot(v(1:2),[0 0],'--r','LineWidth',2)
% subplot(3,1,2)
% errorbar(lon,sec_trend*k,sec_plusminus*k,'o', 'MarkerFace','k')
% v=axis;
% hold on
% plot(v(1:2),[0 0],'--r','LineWidth',2)
% %subplot(4,1,3)
% %errorbar(slon,sec2_trend*k,sec2_plusminus*k,'o', 'MarkerFace','k')
% subplot(3,1,3)
% errorbar(lon,euc_trend*k,euc_plusminus*k,'o', 'MarkerFace','k')
% v=axis;
% hold on
% plot(v(1:2),[0 0],'--r','LineWidth',2)



 
 