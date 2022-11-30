% EUC SODA analysis

close all


% var_ls = {'tas','uas','vas', 'pr','mrso','tos'};
% model_ls = {'ukmo_hadcm3','mpi_echam5','gfdl_cm21', 'pr','mrso','tos'};

%tos_ukmo_hadcm3.nc

%eval (
    ncid = netcdf.open ('tos_csiro_mk35.nc','NC_NOWRITE');
    
    
    
    % get variable ID
    varid_lat = netcdf.inqVarID(ncid,'lat');
    varid_lon = netcdf.inqVarID(ncid,'lon');
    varid_time = netcdf.inqVarID(ncid,'time');
    varid_u = netcdf.inqVarID(ncid,'tos');
    
    
    % get variables
    time = netcdf.getVar(ncid,varid_time);
    lat = netcdf.getVar(ncid,varid_lat);
    lon = netcdf.getVar(ncid,varid_lon);
    u = netcdf.getVar(ncid,varid_u);
   
    
 % lon x lat x depth x time
 
 % eliminate missing value fill value (-9.99e+33)
 % u(u<-10) = NaN;
 
 %% Latitude confinement
 % Indexing spatial bounds 
 I_lat = find(lat<=1 & lat>=-1);
 I_lon = find(lon>180);
 I_dep = find(depth>=25);
 
%  slon = lon(I_lon);
%  slat = lat(I_lat);
%  sdep = depth(I_dep);
 
 % length of various bounds
%  nt = length(time);
%  nlon = length(slon);
%  nlat = length(slat);
%  ndep = length(sdep);
%  
 
 % dimensions: lon x lat x depth x time
 % calculate average zonal velocity over all time and I lat
 avg_u = squeeze(mean(mean(u(:,I_lat,:,:),2),4)); 
 
 % plotting
 figure(100);contour(lon,depth,avg_u','LineWidth',2)
 set(gca,'Ydir','reverse')
 
 %%
 %max_u = squeeze(squeeze(squeeze(max(max(mean(u(:,I_lat,:,:),2))))));
 %max_u = squeeze(squeeze(squeeze(max(max(max(u(I_lon,I_lat,:,:)))))));
 
 slon = lon(I_lon);
 slat = lat(I_lat);
 sdep = depth(I_dep);
 
 nt = length(time);
 nlon = length(slon);
 nlat = length(slat);
 ndep = length(sdep);
 
 % lon x lat x depth x time
%[I J K] = ind2sub([nlon nlat ndep],u(:,:,:,:) == max_u);
 

 Istore = [];
 Jstore = [];
 Kstore = [];
 count = [];
 
 for jj = 1:nt
    
    
    [I J K] = ind2sub([nlon nlat ndep],find(u(I_lon,I_lat,I_dep,jj) == squeeze(squeeze(squeeze(max(max(max(u(:,I_lat,I_dep,jj)))))))));
    %[I J K] = ind2sub([nlon nlat ndep],find(u(:,:,:,jj) == squeeze(squeeze(squeeze(max(max(max(u(:,:,:,jj)))))))));
    
%     I = find(u == squeeze(squeeze(squeeze(max(max(max(u(:,:,:,jj))))))));
    count(jj)=numel(I);
    
    Istore = [Istore; I];
    Jstore = [Jstore; J];
    Kstore = [Kstore; K];
 end
 %%
 
 
 %close all
 figure(2)
 plot3(slon(Istore),slat(Jstore),sdep(Kstore),'o')
 axis([120 280 -5 5 0 300])
 
 set(gca,'ZDir','reverse')
 %%
 

 % determine running meand for a particular time window (in months)
%  win = 1*12+1;    %   
%  win2 = 1*12+1;
%  
%  avgmax_u = runmean(max_u,win);
%  avgmax_u_2 = runmean(max_u,win2);
%  
%  timeser = avgmax_u;
%  
 %Making Trends for particular "modes" in decadal variability 
 %%%%% define trend bounds
%  a = find(time<230,1,'last');
%  b = find(time>285 & time<426.5,1,'first');
%  c = find(time>295 & time<426.5,1,'last');
%  d = find(time>446.5,1,'first');
 
 %%%%%% calculating trends
% trend_1= timeser((1+floor(win/2)):(end-floor(win/2)))-detrend(timeser((1+floor(win/2)):(end-floor(win/2)))); % the trend is the difference between the original time series and the detrended time series
%  trend_2= timeser(b:c)-detrend(timeser(b:c));
%  trend_3 = timeser(d:558)-detrend(timeser(d:558));
%  trend_4 = timeser(43:558)-detrend(timeser(43:558));

%%
 % putting time in datenum format
 %time_2 = datenum(1960,(1+double(time)),1);
 
 time_2 = zeros(1,length(time));
 
 timen = time-.5;
 for j = 1:length(time)
    time_2(j) = addtodate(datenum(1960,1,1),double(timen(j)),'month');
 end
 
 % plot it up
 close all
 
 max_u = squeeze(squeeze(squeeze(max(max(max(u(I_lon,I_lat,I_dep,:)))))));
 max_u_2 = squeeze(squeeze(squeeze(nanmean(max(max(u(I_lon,I_lat,I_dep,:),[],2),[],3)))));
 
 
 titles = {'monthly data','1-yr mean','7-yr mean','10-yr mean'};
 win_i = [0,1,7,10]; % number of years in running mean
 
 close all
  for j = 1:4
 
      % get running mean  
      win = win_i(j)*12+1;
      avgmax_u = runmean(max_u_2,win);timeser = avgmax_u;
      
      % calculate trend
      trend= timeser((1+floor(win/2)):(end-floor(win/2)))-detrend(timeser((1+floor(win/2)):(end-floor(win/2))));
      
     subplot(2,2,j)
     
     % plot data or running mean
     plot(time_2,avgmax_u)
     hold on
     % plot trendline
     plot(time_2((1+floor(win/2)):(end-floor(win/2))),trend,'-.r','LineWidth',2)
     
     %plot details
     legend(char(titles(j)),'Location', 'northwest')
     title((trend(2)-trend(1))*12*100); % this calculates the trend per century if the data are daily
     datetick('x','yyyy','keeplimits')
    
 end
 
 
%  subplot 
%  plot(time_2,avgmax_u_2)
%  hold on
%  %plot(time_2,avgmax_u,'k','LineWidth',2)
%  plot(time_2((1+floor(win/2)):(end-floor(win/2))),trend_1,'-.r','LineWidth',2)
%  %plot(time_2(b:c),trend_2,'k')
%  %plot(time_2(d:558),trend_3,'k')
%  %plot(time_2(43:558),trend_4,'g', 'LineWidth',2)
%  datetick('x','yyyy','keeplimits')
%  ylim([1 1.7])
%  
% title((trend_1(2)-trend_1(1))*12*100); % this calculates the trend per century if the data are daily
% 
% legend('1-yr')


%% Determining EUC Current transport Across -1<=Lat<=1
close all
% calculating area of grid box:
% depth dimension
 I_lat = find(lat<=2 & lat>=-2);

depth_int = diff([0;depth]);
height = zeros(length(depth),1);

for j = 2:length(depth)
    height(j-1) = 0.5*(depth_int(j) + depth_int(j-1));  
end
height(1) = height(1)+ 0.5*depth_int(1);
height(end) = depth_int(end);

% latitudinal dimension: always about 0.5 degrees apart, assume 110.574 km
% per degree lat
% ASSUMPTION: distant between latitudes is constant - not actually true

y_dim = 0.5*110.574*100; %in meters

% calculate grid of areas (square meters)

grid_area_1 = repmat(height*y_dim,1,numel(find(I_lat)));
grid_area_2 = repmat(grid_area_1,[1 1 length(lon)]);
grid_area_3 = permute(grid_area_2,[3 2 1]);
grid_area_4 = repmat(grid_area_3,[1 1 1 length(time)]);

% apply to all SODA velocities between -1<=Lat<=1

transport_1 = grid_area_4.*u(:,(lat<=2 &lat>=-2),:,:);

% only interested in integrating eastward flow (positive values)
% set all westward values (negative) = 0 so only summing east flow
transport_2 = transport_1;
transport_2(transport_1<0) = 0;

% looking at westward flow
transport_3 = transport_1;
transport_3(transport_1>0) = 0;

% looking at net transport in domain -1<=Lat<=1

domain_transport = squeeze(squeeze(squeeze(nansum(nansum(nansum(transport_1))))))/size(transport_1,1);
east_transport = squeeze(squeeze(squeeze(nansum(nansum(nansum(transport_2(:,:,:,:)))))))/size(transport_2,1);
west_transport = squeeze(squeeze(squeeze(nansum(nansum(nansum(transport_3(:,:,(depth>200),:)))))))/size(transport_3,1);
west_transport_2 = -1*west_transport;

lon_domain_transport = squeeze(squeeze(nansum(nansum(transport_1,2),3)))/(10^6);
lon_east_transport = squeeze(squeeze(nansum(nansum(transport_2,2),3)))/(10^6);

%%
% Putting time in datenum form
time_2 = zeros(1,length(time));
 
 timen = time-.5;
 for j = 1:length(time)
    time_2(j) = addtodate(datenum(1960,1,1),double(timen(j)),'month');
 end

titles = {'monthly data','1-yr mean','7-yr mean','10-yr mean'};
 win_i = [0,1,7,10]; % number of years in running mean
 
  for j = 1:4
 
      % get running mean  
      win = win_i(j)*12+1;
      avg_transport = runmean(west_transport_2,win)/(10^6);timeser = avg_transport;
      
      % calculate trend
      trend= timeser((1+floor(win/2)):(end-floor(win/2)))-detrend(timeser((1+floor(win/2)):(end-floor(win/2))));
      
     subplot(2,2,j)
     
     % plot data or running mean
     plot(time_2,avg_transport)
     hold on
     % plot trendline
     plot(time_2((1+floor(win/2)):(end-floor(win/2))),trend,'-.r','LineWidth',2)
     
     %plot details
     legend(char(titles(j)),'Location', 'northwest')
     title((trend(2)-trend(1))*12*100); % this calculates the trend per century if the data are daily
     datetick('x','yyyy','keeplimits')
    
  end
 
%%
 
trend_1 = -500*ones(length(lon),1);
plusminus = trend_1;
sig = trend_1;
trend_2 = trend_1;
plusminus_2 = trend_1;
sig_2 = trend_1;

for j = 1: length(lon)
    [trend_1(j),plusminus(j),sig(j)]=trend_stat(lon_domain_transport(j,:)',99);
    [trend_2(j),plusminus_2(j),sig_2(j)]=trend_stat(lon_east_transport(j,:)',99);  
end

close all
subplot(2,1,1)
errorbar(lon,trend_1,plusminus,'o', 'MarkerFace','k')
subplot(2,1,2)
errorbar(lon,trend_2,plusminus_2,'o', 'MarkerFace','k')

%%
AA = sum(sig)/length(lon);
BB = sum(sig_2)/length(lon);

T1 = trend_1(sig==1);
T2 = trend_2(sig==1);

percent_pos_1 = numel(T1(T1>0))/sum(sig);
percent_pos_2 = numel(T2(T2>0))/sum(sig_2);

%% Seasonal Spread Sheets: Transports

domain_season = double((reshape(domain_transport,12,138)')/(10^6));
eastward_season = double((reshape(east_transport,12,138)')/(10^6));

save domain_season domain_season -ASCII
save eastward_season eastward_season -ASCII

%% Seasonal Spread Sheets: Max Zonal Velocity

maxu_season = double(reshape(max_u,12,138)');
maxu_lon_season = double(reshape(max_u_2,12,138)');

save maxu_season maxu_season -ASCII
save maxu_lon_season maxu_lon_season -ASCII

%% Seasonal Spread Sheets: Average Max Zonal Velocity
