close all 
clear all

 ncid = netcdf.open ('u_eqpac_SODA_extended.nc','NC_NOWRITE'); 
    
    % get variable IDs for SEC
    varid_lat = netcdf.inqVarID(ncid,'lat');
    varid_lon = netcdf.inqVarID(ncid,'lon');
    varid_time = netcdf.inqVarID(ncid,'time');
    varid_depth = netcdf.inqVarID(ncid,'depth');
    
    lat = netcdf.getVar(ncid,varid_lat);
    lon = netcdf.getVar(ncid,varid_lon);
    depth1 = netcdf.getVar(ncid,varid_depth);
    time = netcdf.getVar(ncid,varid_time);
    %%
    lat2 = [(-6.75:.5:-5.25),lat',(5.25:.5:6.75)]';
        

%     ncid = netcdf.open ('u_eqpac_SODA_extended_lat10.nc','NC_NOWRITE');
%     ncid_2 = netcdf.open ('taux_eqpac_SODA_extended_lat10.nc','NC_NOWRITE');    
%     % for looking at EUC too
%     ncid_3 = netcdf.open ('u_eqpac_SODA_extended.nc','NC_NOWRITE'); 
%     
%     % get variable IDs for SEC
%     varid_lat_1 = netcdf.inqVarID(ncid,'lat');
%     varid_lon = netcdf.inqVarID(ncid,'lon');
%     varid_time = netcdf.inqVarID(ncid,'time');
%     varid_u_sec = netcdf.inqVarID(ncid,'u');
%     varid_depth_sec = netcdf.inqVarID(ncid,'depth');
%     
%     % get variable IDs for taux
%     varid_taux = netcdf.inqVarID(ncid_2,'taux');
%     
%     % get variable IDs for EUC
%     varid_u_euc = netcdf.inqVarID(ncid,'u');
%     varid_depth_euc = netcdf.inqVarID(ncid,'depth');
%     varid_lat_2 = netcdf.inqVarID(ncid,'lat');
%     
%     % get variables (general)
%     time = netcdf.getVar(ncid,varid_time);
%     sec_lat = netcdf.getVar(ncid,varid_lat_1);
%     euc_lat = netcdf.getVar(ncid_3,varid_lat_2);
%     lon = netcdf.getVar(ncid,varid_lon);
%     
%     % variables for SEC (lon x lat x depth x time)
%     u_sec = netcdf.getVar(ncid,varid_u_sec);
%     depth_sec = netcdf.getVar(ncid,varid_depth_sec);
%     
%     % variables for taux (lon x lat x time)
%     taux = netcdf.getVar(ncid_2,varid_taux);
%     
%     % variables for EUC (lon x lat x depth x time)
%     
%     u_euc = netcdf.getVar(ncid_3,varid_u_euc);
%     depth_euc = netcdf.getVar(ncid_3,varid_depth_euc);
%  
%  % eliminate missing value fill value (-9.99e+33)
%  u_sec(u_sec<-10000) = NaN;
%  taux(taux<-10000) = NaN;
%  u_euc(u_euc<-10000) = NaN;
%  

 %% Indexing spatial bounds
 
 I_lat_1 = find(lat2<=7 & lat2>=-7);      % For SEC & taux
 I_lat_2 = find(lat2<=2 & lat2>=-2);      % For EUC
 I_lon = find(lon>=160 &lon<=260);
 
 slat_1 = lat2(I_lat_1);
 slat_2 = lat2(I_lat_2);
 slon = lon(I_lon);

 nt = length(time);
 nlon = length(slon);
 nlat_1 = length(slat_1);
 nlat_2 = length(slat_2);
  
 %%

 load EUC_data.mat
 load U_EUC.mat
 load U_SEC.mat
 load TAUX.mat
 
 %%
 
% U_EUC = squeeze(mean(u_euc(I_lon,I_lat_2,:,:),2));
% U_SEC = squeeze(u_sec(I_lon,I_lat_1,:,:)); 
% TAUX = squeeze(taux(I_lon,I_lat_1,:));
 
 
%% Climatology for winds and currents...

EUCp = permute(U_EUC,[3,2,1]);
[EUC_clim EUC_anom] =  climanom(EUCp);
SECp = permute(U_SEC,[3,2,1]);
[SEC_clim SEC_anom] =  climanom(SECp);
TAUXp = permute(TAUX,[3,2,1]);
[TAUX_clim TAUX_anom] =  climanom(TAUXp);
% 
% %%
% TEMPp = permute(TEMP,[3,2,1]);
% [TEMP_clim TEMP_anom] =  climanom(TEMPp);
% 
% SLPp = permute(double(SLP),[3,2,1]);
% [SLP_clim SLP_anom] =  climanom(SLPp);
% 
% SSHp = permute(SSH,[3,2,1]);
% [SSH_clim SSH_anom] =  climanom(SSHp);
% 
% SSTp = permute(SST,[3,2,1]);
% [SST_clim SST_anom] =  climanom(SSTp);
%%
% U_SEC_2 = permute(SEC_anom,[3,2,1]);
% TAUX_2 = permute(TAUX_anom,[3,2,1]);
% U_EUC_2 = permute(EUC_anom,[3,2,1]);
% 
% %% Trends in anomaly
% taux_trend = zeros(12,length(slon),length(lat));
% taux_plusminus = taux_trend;
% taux_sig = taux_trend;
% 
% sec_trend = zeros(12,length(slon),length(lat));
% sec_plusminus = sec_trend;
% sec_sig = sec_trend;
% 
% euc_trend = zeros(12,length(slon),length(depth));
% euc_plusminus = euc_trend;
% euc_sig = euc_trend;
% 
% % for every month 
% for j = 1:12
%     nmon = 12*(1:137)+j;
%     % for every lon
%     for k = 1: length(slon)
%         
%         % for every latitude
%         for m = 1:length(lat)
%             [sec_trend(j,k,m),sec_plusminus(j,k,m),sec_sig(j,k,m)]=trend_stat(squeeze(U_SEC_2(k,m,nmon)),99);    
%             [taux_trend(j,k,m),taux_plusminus(j,k,m),taux_sig(j,k,m)]=trend_stat(squeeze(TAUX_2(k,m,nmon)),99);
%         end
%         
%         % for every depth
%         for n = 1:length(depth)
%             [euc_trend(j,k,n),euc_plusminus(j,k,n),euc_sig(j,k,n)]=trend_stat(squeeze(U_EUC_2(k,n,nmon)),99);
%         end
%     end
% end

load TREND_SEC.mat
load TREND_EUC.mat
load TREND_TAUX.mat
% load TREND_TEMP.mat
% load TREND_SLP.mat
% load TREND_SSH.mat
% load TREND_SST.mat
%%
close all

% set figure window size
fig1 = figure('Color',[1 1 1]);
ss=get(0,'ScreenSize');
set(fig1,'Position',[ss(1)+15 ss(2)+25 ss(3)-200 ss(4)-80]) 
winsize = get(fig1,'Position');
winsize(1:2) = [0 0];
numframes=12;
A=moviein(numframes,fig1,winsize);
set(fig1,'NextPlot','replacechildren')
month = {'January', 'February','March','April','May','June','July','August','September','October','November','December'}; 

%mTauxT = taux_trend - nanmean(nanmean(nanmean(taux_trend)));
%mSecT = sec_trend - nanmean(nanmean(nanmean(sec_trend)));
%mEucT = euc_trend - nanmean(nanmean(nanmean(euc_trend)));

% mTauxC = TAUX_clim - nanmean(nanmean(nanmean(TAUX_clim)));
% mSecC = SEC_clim - nanmean(nanmean(nanmean(SEC_clim)));
% mEucC = EUC_clim - nanmean(nanmean(nanmean(EUC_clim)));

% mTempT = temp_trend - nanmean(nanmean(nanmean(temp_trend)));
% mSstT = sst_trend - nanmean(nanmean(nanmean(sst_trend)));
% mSlpT = slp_trend - nanmean(nanmean(nanmean(slp_trend)));
% mSshT = ssh_trend - nanmean(nanmean(nanmean(ssh_trend)));

% mTempC = TEMP_clim - nanmean(nanmean(nanmean(TEMP_clim)));
% mSstC = SST_clim - nanmean(nanmean(nanmean(SST_clim)));
% mSlpC = SLP_clim - nanmean(nanmean(nanmean(SLP_clim)));
% mSshC = SSH_clim - nanmean(nanmean(nanmean(SSH_clim)));

for j = 1
   
clf
 
% get the corners of the domain in which the data occurs.
min_x = min(min(slon));
min_y = min(min(lat2));
max_x = max(max(slon));
max_y = max(max(lat2));

min_y3 = min(min(lat));
max_y3 = max(max(lat));

min_y4= min(min(slat_3));
max_y4 = max(max(slat_3));

min_y2 = max_y + 10;
min_y2 = min_y2 +(max_y-min_y);
 
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
imgzposition = 400;         % Winds 
imgzposition2 = 180;        % SEC
imgzposition3 = 0;          % EUC

%%%%%%%%%%%%%%%%%%%%  ADDING TEMP, SST, SLP, SSH  %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % normalized trends here
% planeimg2_1 = squeeze(temp_trend(j,:,:))/(max(max(max(abs(temp_trend)))));
% planeimg2_2 = squeeze(sst_trend(j,:,:))/(max(max(max(abs(sst_trend)))));
% planeimg2_3 = squeeze(slp_trend(j,:,:))/(max(max(max(abs(slp_trend))))); 
% planeimg2_4 = squeeze(ssh_trend(j,:,:))/(max(max(max(abs(ssh_trend))))); 
% 
% % normalized climatology values
% 
% planeimg2_5 = squeeze(mTempC(j,:,:))/(max(max(max(abs(mTempC)))));
% planeimg2_6 = squeeze(mSstC(j,:,:))/(max(max(max(abs(mSstC)))));
% planeimg2_7 = squeeze(mSlpC(j,:,:))/(max(max(max(abs(mSlpC))))); 
% planeimg2_8 = squeeze(mSshC(j,:,:))/(max(max(max(abs(mSshC))))); 
% 
% % desired z position of the image plane.
% imgzposition2_1 = 0;      % TEMP
% imgzposition2_2 = 120;      % SST
% imgzposition2_3 = 350;        % SLP
% imgzposition2_4 = 350;        % SSH  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% plot plane that represents the ocean surface. code it to be blue (i.e. *-2)
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
%[C h] = contourf(slon,lat2,planeimg1',100,'LineWidth',2);
[C h] = contourf(slon,lat2,planeimg4,100,'LineWidth',2);
set(h, 'EdgeColor','none')

% TAUX climatology contours
%[C p] = contour(slon,lat2,planeimg4,5,'LineWidth',2);


%rotate plots into appropriate dimension
zdir = [0 1 0];
origin1 = [0 0 imgzposition];
rotate(get(h,'children'),zdir,0,origin1) 
%rotate(get(p,'children'),zdir,0,origin1)


%%% TRADES plot (first layer-shift)
% dark background surface for contrast
% surf([max_x (max_x+(max_x-min_x))]+.5,[max_y (max_y+(max_y-min_y))]+.5,repmat(-1*imgzposition+.5, [2 2]),...
%        NaN(size(planeimg4)),'facecolor','texture','FaceAlpha',1)
% 
% % TAUX trend filled contours  
% [C h] = contourf(slon,lat2,planeimg1',100,'LineWidth',2);
% set(h, 'EdgeColor','none')
% 
% % TAUX climatology contours
% [C p] = contour(slon,lat2,planeimg4,5,'LineWidth',2);
% 
% %rotate plots into appropriate dimension
% zdir = [0 1 0];
% origin1 = [0 0 imgzposition];
% 
% nxdat = get(h,'Xdata')+(max_x-min_x)+.5;
% set(h,'Xdata',nxdat)
% nydat = get(h,'Ydata')+(max_y-min_y)+.5;
% set(h,'Ydata',nydat)
% 
% rotate(get(h,'children'),zdir,0,origin1) 
% 
% 
% nxdat = get(p,'Xdata')+(max_x-min_x)+.5;
% set(p,'Xdata',nxdat)
% nydat = get(p,'Ydata')+(max_y-min_y)+.5;
% set(p,'Ydata',nydat)
% rotate(get(p,'children'),zdir,0,origin1)

% %%% SSH plot (second layer)
% % dark background surface for contrast
% %surf([min_x max_x],[min_y3 max_y3],repmat(-1*imgzposition2_4+.5, [2 2]),...
%        NaN(size(planeimg2_8)),'facecolor','texture','FaceAlpha',1)
% 
% % SSH trend filled contours  
% [C h] = contourf(slon,lat,planeimg2_4',100,'LineWidth',2);
% set(h, 'EdgeColor','none')
% 
% % SSH climatology contours
% [C p] = contour(slon,lat,planeimg2_8,7,'LineWidth',2);
% 
% %rotate plots into appropriate dimension
% zdir = [0 1 0];
% origin1 = [0 0 imgzposition2_4];
% rotate(get(h,'children'),zdir,0,origin1) 
% rotate(get(p,'children'),zdir,0,origin1)




% %%% SLP plot (second layer-shift)
% % dark background surface for contrast
% surf([max_x (max_x+(max_x-min_x))]+.5,[min_y4 max_y4]+(max_y-min_y+0.5),repmat(-1*imgzposition2_4+.5, [2 2]),...
%        NaN(size(planeimg2_8)),'facecolor','texture','FaceAlpha',1)
% 
% % SLP trend filled contours  
% [C h] = contourf(SLP_lon,slat_3,planeimg2_3',100,'LineWidth',2);
% set(h, 'EdgeColor','none')
% 
% % SLP climatology contours
% [C p] = contour(SLP_lon,slat_3,planeimg2_7,7,'LineWidth',2);
% 
% %rotate plots into appropriate dimension
% zdir = [0 1 0];
% origin1 = [0 0 imgzposition2_4];
% 
% nxdat = get(h,'Xdata')+(max_x-min_x)+.5;
% set(h,'Xdata',nxdat)
% nydat = get(h,'Ydata')+(max_y-min_y)+.5;
% set(h,'Ydata',nydat)
% rotate(get(h,'children'),zdir,0,origin1) 
% 
% nxdat = get(p,'Xdata')+(max_x-min_x)+.5;
% set(p,'Xdata',nxdat)
% nydat = get(p,'Ydata')+(max_y-min_y)+.5;
% set(p,'Ydata',nydat)
% rotate(get(p,'children'),zdir,0,origin1)





%%% SEC plot (third layer)
% dark background surface for contrast
surf([min_x max_x],[min_y max_y],repmat(-1*imgzposition2+.1, [2 2]),...
      NaN(size(planeimg4)),'facecolor','texture','FaceAlpha',1)
 
% SEC trend filled contours
%[C h] = contourf(slon,lat2,planeimg2',100,'LineWidth',2);
[C h] = contourf(slon,lat2,planeimg5,100,'LineWidth',2);
set(h, 'EdgeColor','none')

% SEC climatology contours
%[C p] = contour(slon,lat2,planeimg5,5,'LineWidth',2);

%rotate plots into appropriate dimension
zdir = [0 1 0];
origin1 = [0 0 imgzposition2];
rotate(get(h,'children'),zdir,0,origin1) 
%rotate(get(p,'children'),zdir,0,origin1)

% 
% 
% %%% SST plot (third layer - shift)
% % dark background surface for contrast
%  surf([max_x (max_x+(max_x-min_x))]+.5,([min_y3 max_y3] +(max_y-min_y+0.5)),repmat(-1*imgzposition2+.1, [2 2]),...
%        NaN(size(planeimg2_6)),'facecolor','texture','FaceAlpha',1);
%  
% % SST trend filled contours
% [C h] = contourf(slon,lat,planeimg2_2',100,'LineWidth',2);
% set(h, 'EdgeColor','none')
% 
% % SST climatology contours
% [C p] = contour(slon,lat,planeimg2_6,5,'LineWidth',2);
% 
% %rotate plots into appropriate dimension
% zdir = [0 1 0];
% origin1 = [0 0 imgzposition2];
% % origin2 = [12.5 12.5 imgzposition2];
% 
% 
% nxdat = get(h,'Xdata')+(max_x-min_x)+.5;
% set(h,'Xdata',nxdat)
% nydat = get(h,'Ydata')+(max_y-min_y)+.5;
% set(h,'Ydata',nydat)
% 
% rotate(get(h,'children'),zdir,0,origin1) 
% 
% 
% nxdat = get(p,'Xdata')+(max_x-min_x)+.5;
% set(p,'Xdata',nxdat)
% nydat = get(p,'Ydata')+(max_y-min_y)+.5;
% set(p,'Ydata',nydat)
% rotate(get(p,'children'),zdir,0,origin1)




%%% EUC PLOT (bottom layer)

h12 = surf([min_x max_x],[min(depth1) max(depth1)],repmat(-1*imgzposition3-.1, [2 2]),...
       NaN(size(planeimg6)),'facecolor','texture','FaceAlpha',1);

% EUC trend filled contours
%[C h] = contourf(slon,depth1,planeimg3',100,'LineWidth',2);
[C h] = contourf(slon,depth1,planeimg6,100,'LineWidth',2);
set(h, 'EdgeColor','none')

% EUC climatology contours
%[C p] = contour(slon,depth1,planeimg6,5,'LineWidth',2);

%rotate plots into appropriate dimension
zdir = [1 0 0];
origin1 = [0 0 imgzposition3];
rotate(h12,zdir,90,origin1)
rotate(get(h,'children'),zdir,90,origin1) 
%rotate(get(p,'children'),zdir,90,origin1)
set(gca,'ZDir','reverse')


% 
% %%% TEMP PLOT (bottom layer - shift)
% 
% h12 = surf([min_x max_x]+(max_x-min_x+2),[min(depth1) max(depth1)],repmat(-1*imgzposition3-.1, [2 2]),...
%        NaN(size(planeimg6)),'facecolor','texture','FaceAlpha',1);
% 
% % TEMP trend filled contours
% [C h] = contourf(slon,depth1,planeimg2_1',100,'LineWidth',2);
% set(h, 'EdgeColor','none')
% 
% % TEMP climatology contours
% [C p] = contour(slon,depth1,planeimg2_5,9,'LineWidth',2);
% 
% %shift and rotate plots into appropriate dimension
% zdir = [1 0 0];
% origin1 = [0 0 imgzposition3];
% origin2 = [12.5 12.5 (imgzposition3)];
% rotate(h12,zdir,90,origin2)
% 
% nxdat = get(h,'Xdata')+(max_x-min_x+2);
% set(h,'Xdata',nxdat)
% rotate(get(h,'children'),zdir,90,origin2)
% 
% nxdat = get(p,'Xdata')+(max_x-min_x+2);
% set(p,'Xdata',nxdat)
% rotate(get(p,'children'),zdir,90,origin2)
% 
% set(gca,'ZDir','reverse')

%%% Finishing Touches %%%%%%%%%%%%%

% plotting dashed lines (guides for spatial orientation)
 plot3([max_x,max_x],[0,0],[min(depth1), -1*imgzposition],'--k','linewidth',2)
 plot3([min_x,min_x],[0,0],[min(depth1), -1*imgzposition],'--k','linewidth',2)
 
 %plot3([max_x,max_x],[(max_y-min_y+0.5),(max_y-min_y+0.5)],[min(depth1), -1*imgzposition],'--k','linewidth',2)
 %plot3([max_x+(max_x-min_x),max_x+(max_x-min_x)]+.5,[(max_y-min_y+0.5),(max_y-min_y+0.5)],[min(depth1), -1*imgzposition],'--k','linewidth',2)

 % plotting equators (guides for spatial orientation)
 plot3([min_x,max_x],[0,0],[-1*imgzposition+.1,-1*imgzposition+.1],'-k','linewidth',2)
%  plot3([min_x,max_x],[0,0],[-1*imgzposition2_4+.1,-1*imgzposition2_4+.1],'-k','linewidth',2)
 plot3([min_x,max_x],[0,0],[-1*imgzposition2+.1,-1*imgzposition2+.1],'-k','linewidth',2)
 plot3([min_x,max_x],[0,0],[min(depth1), min(depth1)],'-k','linewidth',2)
 plot3([min_x,max_x],[0,0],[max(depth1), max(depth1)],'-k','linewidth',2)

%  plot3([min_x,max_x]+(max_x-min_x+0.5),[(max_y-min_y+0.5),(max_y-min_y+0.5)],[-1*imgzposition+.1,-1*imgzposition+.1],'-k','linewidth',2)
%  % plot3([min_x,max_x]+(max_x-min_x+0.5),[(max_y-min_y+0.5),(max_y-min_y+0.5)],[-1*imgzposition2_4+.1,-1*imgzposition2_4+.1],'-k','linewidth',2)
%  plot3([min_x,max_x]+(max_x-min_x+0.5),[(max_y-min_y+0.5),(max_y-min_y+0.5)],[-1*imgzposition2+.1,-1*imgzposition2+.1],'-k','linewidth',2)
%  plot3([min_x,max_x]+(max_x-min_x+0.5),[(max_y-min_y+0.5),(max_y-min_y+0.5)],[min(depth1), min(depth1)],'-k','linewidth',2)
%  plot3([min_x,max_x]+(max_x-min_x+0.5),[(max_y-min_y+0.5),(max_y-min_y+0.5)],[max(depth1), max(depth1)],'-k','linewidth',2)
 
% plotting equator lines at depth (guides spatial orientation)
 plot3([max_x,max_x],[0,0],[min(depth1), max(depth1)],'-k','linewidth',2)
% plot3([min_x,min_x],[0,0],[min(depth1), max(depth1)],'-k','linewidth',2)

% set axis limits
xlim([160 260])

zlim([-450 max(depth1)])
 
% set the view angle.
view(20,15);
whitebg([.8 .8 .8])
% colorbar('YTickLabel')
colorbar('YTickLabel',[],'YTick',[])

% LABELS
%title(char(month(j)),'FontSize',30,'Color',[0 0 0])

% xlabel('Longitude','FontSize',25);
% ylabel('Latitude','FontSize',25);
% zlabel('Depth (m)','FontSize',25);

 z=0:100:600;
 set(gca,'ZTick',z)
 %set(gca,'ZTickLabel',[' 0 ';'100';'200';'300';'400';'500';'600']);
 set(gca,'ZTickLabel',[]);
 set(gca,'fontsize',20)
 
 x = 160:20:260;
 set(gca,'XTick',x)
 %set(gca,'XTickLabel',['160E';'180 ';'160W';'140W';'120W';'100W']);
 set(gca,'XTickLabel',[]);
 set(gca,'fontsize',20)
 
 y = -4:4:4;
 set(gca,'YTick',y)
 %set(gca,'YTickLabel',['-4 ';' 0 ';' 4 ']);
 set(gca,'YTickLabel',[]);
 set(gca,'fontsize',20)
 
% reset background color to light grey

% custom colormap for the figure
colormap([0.1 0.1 0.1;jet(250)]);
caxis([-1.01 1])

% capture fraome for movie before going on to next month
% A(:,j)=getframe(fig1,winsize);



end

% concatenate video structure to repreat video loop
% A2 = [A,A,A,A,A,A];
% %%
% % watch the movie to make sure it works
% movie(fig1,A,3,1,winsize)
% %save EUC_clim.mat A2 
% 
% %% convert to avi format which can go in power point presentation
% 
% %movie2avi(A2, 'EUC_clim.avi','fps',1,'compression', 'none');
% movie2avi(A2, 'EUC_trends2.avi','fps',1,'compression', 'none');
