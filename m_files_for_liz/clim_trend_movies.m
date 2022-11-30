close all 
clear all

tic
 ncid = netcdf.open ('u_eqpac_SODA_extended.nc','NC_NOWRITE'); 
    
    % get variable IDs for SEC
    varid_lat = netcdf.inqVarID(ncid,'lat');
    varid_lon = netcdf.inqVarID(ncid,'lon');
    varid_time = netcdf.inqVarID(ncid,'time');
    varid_depth = netcdf.inqVarID(ncid,'depth');
    
    lat = netcdf.getVar(ncid,varid_lat);
    full_lon = netcdf.getVar(ncid,varid_lon);
    depth = netcdf.getVar(ncid,varid_depth);
    time = netcdf.getVar(ncid,varid_time);
    
    lat = [(-6.75:.5:-5.25),lat',(5.25:.5:6.75)]';
        

    ncid = netcdf.open ('u_eqpac_SODA_extended_lat10.nc','NC_NOWRITE');
    ncid_2 = netcdf.open ('taux_eqpac_SODA_extended_lat10.nc','NC_NOWRITE');    
    % for looking at EUC too
    ncid_3 = netcdf.open ('u_eqpac_SODA_extended.nc','NC_NOWRITE'); 
    
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
 I_lat_2 = find(euc_lat<=3 & euc_lat>=-3);      % For EUC
 I_lon = find(lon>=160 &lon<=260);
 
 slat_1 = sec_lat(I_lat_1);
 slat_2 = euc_lat(I_lat_2);
 slon = lon(I_lon);

 nt = length(time);
 nlon = length(slon);
 nlat_1 = length(slat_1);
 nlat_2 = length(slat_2);
  
 %%
 
U_EUC = squeeze(mean(u_euc(I_lon,I_lat_2,:,:),2));
U_SEC = squeeze(u_sec(I_lon,I_lat_1,:,:)); 
TAUX = squeeze(taux(I_lon,I_lat_1,:));
 
 
%% Climatology for winds and currents...

EUCp = permute(U_EUC,[3,2,1]);
[EUC_clim EUC_anom] =  climanom(EUCp);
SECp = permute(U_SEC,[3,2,1]);
[SEC_clim SEC_anom] =  climanom(SECp);
TAUXp = permute(TAUX,[3,2,1]);
[TAUX_clim TAUX_anom] =  climanom(TAUXp);

U_SEC = permute(SEC_anom,[3,2,1]);
TAUX = permute(TAUX_anom,[3,2,1]);
U_EUC = permute(EUC_anom,[3,2,1]);

%% Trends in anomaly
taux_trend = zeros(12,length(slon),length(lat));
taux_plusminus = taux_trend;
taux_sig = taux_trend;

sec_trend = zeros(12,length(slon),length(lat));
sec_plusminus = sec_trend;
sec_sig = sec_trend;

euc_trend = zeros(12,length(slon),length(depth));
euc_plusminus = euc_trend;
euc_sig = euc_trend;

% for every month 
for j = 1:12
    nmon = 12*(1:137)+j;
    % for every lon
    for k = 1: length(slon)
        
        % for every latitude
        for m = 1:length(lat)
            [sec_trend(j,k,m),sec_plusminus(j,k,m),sec_sig(j,k,m)]=trend_stat(squeeze(U_SEC(k,m,nmon)),99);    
            [taux_trend(j,k,m),taux_plusminus(j,k,m),taux_sig(j,k,m)]=trend_stat(squeeze(TAUX(k,m,nmon)),99);
        end
        
        % for every depth
        for n = 1:length(depth)
            [euc_trend(j,k,n),euc_plusminus(j,k,n),euc_sig(j,k,n)]=trend_stat(squeeze(U_EUC(k,n,nmon)),99);
        end
    end
end

%%
close all

% set figure window size
fig1 = figure('Color',[1 1 1]);
ss=get(0,'ScreenSize');
set(fig1,'Position',[ss(1)+5 ss(2) ss(3)-250 ss(4)-120]) 
winsize = get(fig1,'Position');
winsize(1:2) = [0 0];
numframes=12;
A=moviein(numframes,fig1,winsize);
set(fig1,'NextPlot','replacechildren')
month = {'January', 'February','March','April','May','June','July','August','September','October','November','December'};

for j = 1:12
clf
 
% get the corners of the domain in which the data occurs.
min_x = min(min(slon));
min_y = min(min(lat));
max_x = max(max(slon));
max_y = max(max(lat));
 
% the image data you want to show as a plane (the normalized climatology
% data)

% normalized trends here
planeimg1 = squeeze(taux_trend(j,:,:))/(max(max(max(abs(taux_trend)))));
planeimg2 = squeeze(sec_trend(j,:,:))/(max(max(max(abs(sec_trend)))));
planeimg3 = squeeze(euc_trend(j,:,:))/(max(max(max(abs(euc_trend))))); 

% normalized climatology values
planeimg4 = squeeze(TAUX_clim(j,:,:))/(max(max(max(abs(TAUX_clim)))));
planeimg5 = squeeze(SEC_clim(j,:,:))/(max(max(max(abs(SEC_clim)))));
planeimg6 = squeeze(EUC_clim(j,:,:))/(max(max(max(abs(EUC_clim))))); 

% desired z position of the image plane.
imgzposition = 350;
imgzposition2 = 150;
imgzposition3 = 0;

% plot plane that represents the ocean. code it to be blue (i.e. *-2)
 surf([min_x max_x],[min_y max_y],repmat(4, [2 2]),...
     -.2*ones(size(planeimg4)),'facecolor','texture','FaceAlpha',0.1) 

hold on;

% set inital background color to code NaN values in contour plots
whitebg([0.1 0.1 0.1])  % dark grey background


%%% TRADES plot (first layer)
% dark background surface for contrast
surf([min_x max_x],[min_y max_y],repmat(-1*imgzposition+.5, [2 2]),...
       NaN(size(planeimg4)),'facecolor','texture','FaceAlpha',1)

% TAUX trend filled contours  
[C h] = contourf(slon,lat,planeimg1',100,'LineWidth',2);
set(h, 'EdgeColor','none')

% TAUX climatology contours
[C p] = contour(slon,lat,planeimg4,5,'LineWidth',2);

%rotate plots into appropriate dimension
zdir = [0 1 0];
origin1 = [0 0 imgzposition];
rotate(get(h,'children'),zdir,0,origin1) 
rotate(get(p,'children'),zdir,0,origin1)


%%% SEC plot (second layer)
% dark background surface for contrast
surf([min_x max_x],[min_y max_y],repmat(-1*imgzposition2+.1, [2 2]),...
      NaN(size(planeimg4)),'facecolor','texture','FaceAlpha',1)
 
% SEC trend filled contours
[C h] = contourf(slon,lat,planeimg2',100,'LineWidth',2);
set(h, 'EdgeColor','none')

% SEC climatology contours
[C p] = contour(slon,lat,planeimg5,5,'LineWidth',2);

%rotate plots into appropriate dimension
zdir = [0 1 0];
origin1 = [0 0 imgzposition2];
rotate(get(h,'children'),zdir,0,origin1) 
rotate(get(p,'children'),zdir,0,origin1)


%%% EUC PLOT (bottom layer)

h12 = surf([min_x max_x],[min(depth) max(depth)],repmat(-1*imgzposition3-.1, [2 2]),...
       NaN(size(planeimg6)),'facecolor','texture','FaceAlpha',1);

% EUC trend filled contours
[C h] = contourf(slon,depth,planeimg3',100,'LineWidth',2);
set(h, 'EdgeColor','none')

% EUC climatology contours
[C p] = contour(slon,depth,planeimg6,5,'LineWidth',2);

%rotate plots into appropriate dimension
zdir = [1 0 0];
origin1 = [0 0 imgzposition3];
rotate(h12,zdir,90,origin1)
rotate(get(h,'children'),zdir,90,origin1) 
rotate(get(p,'children'),zdir,90,origin1)
set(gca,'ZDir','reverse')

%%% Finishing Touches %%%%%%%%%%%%%

% plotting dashed lines (guides spatial orientation)
 plot3([max_x,max_x],[0,0],[min(depth), -1*imgzposition],'--k','linewidth',2)
 plot3([min_x,min_x],[0,0],[min(depth), -1*imgzposition],'--k','linewidth',2)

 % plotting equators (guides spatial orientation)
 plot3([min_x,max_x],[0,0],[-1*imgzposition2+.1,-1*imgzposition2+.1],'-k','linewidth',2)
 plot3([min_x,max_x],[0,0],[-1*imgzposition+.1,-1*imgzposition+.1],'-k','linewidth',2)
 plot3([min_x,max_x],[0,0],[min(depth), min(depth)],'-k','linewidth',2)
 plot3([min_x,max_x],[0,0],[max(depth), max(depth)],'-k','linewidth',2)

% plotting equator lines at depth (guides spatial orientation)
 plot3([max_x,max_x],[0,0],[min(depth), max(depth)],'-k','linewidth',2)
 plot3([min_x,min_x],[0,0],[min(depth), max(depth)],'-k','linewidth',2)

% set axis limits
caxis([-1 1]) 
xlim([160 260])
zlim([-400 max(depth)])
 
% set the view angle.
view(20,14);

% LABELS
title(char(month(j)),'FontSize',30,'Color',[0 0 0 ])
% xlabel('Longitude','FontSize',25);
% ylabel('Latitude','FontSize',25);
% zlabel('Depth (m)','FontSize',25);

 z=0:100:600;
 set(gca,'ZTick',z)
 set(gca,'ZTickLabel',[' 0 ';'100';'200';'300';'400';'500';'600']);
 set(gca,'fontsize',20)
 
% reset background color to light grey
whitebg([.8 .8 .8])

% custom colormap for the figure
colormap([0.1 0.1 0.1;jet]);
colorbar

% capture fraome for movie before going on to next month
A(:,j)=getframe(fig1,winsize);

end

% concatenate video structure to repreat video loop
A2 = [A,A,A,A,A,A];

% watch the movie to make sure it works
movie(fig1,A,3,1,winsize)
%save EUC_clim.mat A2 

%% convert to avi format which can go in power point presentation

%movie2avi(A2, 'EUC_clim.avi','fps',1,'compression', 'none');
movie2avi(A2, 'EUC_trends.avi','fps',1,'compression', 'none');
toc
