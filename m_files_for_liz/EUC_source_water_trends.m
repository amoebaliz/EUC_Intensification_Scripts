% close all 
% clear all
% 
% % access data files
% ncid = netcdf.open ('u_eqpac_SODA_extended.nc','NC_NOWRITE'); 
% ncid2 = netcdf.open ('v_eqpac_SODA_extended.nc','NC_NOWRITE');
% ncid3 = netcdf.open ('w_eqpac_SODA_extended.nc','NC_NOWRITE'); 
%     
% % get domain variables 
% % NOTE: assuming lat, lon, depth, time) are the same for u, v, and w 
%     
%     varid_lat = netcdf.inqVarID(ncid,'lat');
%     varid_lon = netcdf.inqVarID(ncid,'lon');
%     varid_time = netcdf.inqVarID(ncid,'time');
%     varid_depth = netcdf.inqVarID(ncid,'depth');
%     
%     lat = netcdf.getVar(ncid,varid_lat);
%     lon = netcdf.getVar(ncid,varid_lon);
%     depth = netcdf.getVar(ncid,varid_depth);
%     time = netcdf.getVar(ncid,varid_time);
%     
%     clear varid_lat varid_lon varid_lat varid_time varid_time
%     
%     % get data variables (u, v, w); Dimensions: <lon x lat x depth x time>
%     varid_u = netcdf.inqVarID(ncid,'u');
%     varid_v = netcdf.inqVarID(ncid2,'v');
%     varid_w = netcdf.inqVarID(ncid3,'w');
%     
%     u = netcdf.getVar(ncid,varid_u);
%     v = netcdf.getVar(ncid2,varid_v);
%     w = netcdf.getVar(ncid3,varid_w);
%     
%     clear varid_u varid_v varid_w ncid ncid2 ncid3 
%         
%     % replace "missing" / fill value (-9.99e+33) with NaN
%     u(u< -10000) = NaN;
%     v(v< -10000) = NaN;
%     w(w< -10000) = NaN;
%     
%  %% Compute Velocity Field Climatology 
%  % NOTE: climanom input dimensions should be <time x lat x lon> so permute
%  % data values, loping through each depth
%  
%  % Initialize u/v/w_clim variables
%  
%  u_clim = zeros(12, length(lat), length(lon), length(depth));
%  v_clim = zeros(12, length(lat), length(lon), length(depth));
%  w_clim = zeros(12, length(lat), length(lon), length(depth));
%  
% % Initialize u/v/w_anom variables
%   
%  u_anom = zeros(length(time), length(lat), length(lon), length(depth));
%  v_anom = zeros(length(time), length(lat), length(lon), length(depth));
%  w_anom = zeros(length(time), length(lat), length(lon), length(depth));
%   
%  
%  for j = 1: length(depth)
%  
%     u_permute = permute(squeeze(u(:,:,j,:)),[3,2,1]);
%     [u_clim(:,:,:,j) u_anom(:,:,:,j)] =  climanom(u_permute);
% 
%     v_permute = permute(squeeze(v(:,:,j,:)),[3,2,1]);
%     [v_clim(:,:,:,j) v_anom(:,:,:,j)] =  climanom(v_permute);
% 
%     w_permute = permute(squeeze(w(:,:,j,:)),[3,2,1]);
%     [w_clim(:,:,:,j) w_anom(:,:,:,j)] =  climanom(w_permute);
% 
%  end
%  
%  clear u_permute v_permute w_permute
%  
%  % NOTE: these fields have dimensions: <time x lat x lon x depth>
%  
%   
%  %% Compute Seasonal (Monthly) Trends in Velocity   
%  % NOTE: Should I be using the anomaly
%  
% % Initialize variables  
% u_clim_trend = zeros(length(lon),length(lat),length(depth),12);
% u_clim_plusminus = u_clim_trend;
% u_clim_sig = u_clim_trend;
% 
% v_clim_trend = zeros(length(lon),length(lat),length(depth),12);
% v_clim_plusminus = v_clim_trend;
% v_clim_sig = v_clim_trend;
% 
% w_clim_trend = zeros(length(lon),length(lat),length(depth),12);
% w_clim_plusminus = w_clim_trend;
% w_clim_sig = w_clim_trend;
%  
% % for every month 
% for mon = 1:12
%      nmon = (12*(0:137))+mon;
%      
%      % for every depth
%      for d = 1: length(depth)
%          
%         % for every longitude
%          for ln = 1: length(lon)
%              
%              % for every latitude
%              for lt = 1:length(lat)
%                  
%                   [u_clim_trend(ln,lt,d,mon) u_clim_plusminus(ln,lt,d,mon), u_clim_sig(ln,lt,d,mon)] = trend_stat(squeeze(u(ln,lt,d,nmon)),99);
%                   [v_clim_trend(ln,lt,d,mon) v_clim_plusminus(ln,lt,d,mon), v_clim_sig(ln,lt,d,mon)] = trend_stat(squeeze(v(ln,lt,d,nmon)),99);
%                   [w_clim_trend(ln,lt,d,mon) w_clim_plusminus(ln,lt,d,mon), w_clim_sig(ln,lt,d,mon)] = trend_stat(squeeze(w(ln,lt,d,nmon)),99);
%     
%              end
%          end
%      end
% end
% 
% % NaN out regions where trends are not significant
% u_clim_trend(u_clim_sig == 0) = NaN;
% v_clim_trend(v_clim_sig == 0) = NaN;
% w_clim_trend(w_clim_sig == 0) = NaN;
% 
%  
%  %% Compute Long-Term (All SODA) Trends in Velocity
% 
%  u_trend = zeros(length(lon),length(lat),length(depth));
%  u_plusminus = u_trend;
%  u_sig = u_trend;
%  
%  v_trend = zeros(length(lon),length(lat),length(depth));
%  v_plusminus = v_trend;
%  v_sig = v_trend;
%  
%  w_trend = zeros(length(lon),length(lat),length(depth));
%  w_plusminus = w_trend;
%  w_sig = w_trend;
%  
%  % for every depth
%  for d = 1: length(depth)
% 
%     % for every longitude
%      for ln = 1: length(lon)
% 
%          % for every latitude
%          for lt = 1:length(lat)
% 
%               [u_trend(ln,lt,d) u_plusminus(ln,lt,d), u_sig(ln,lt,d)] = trend_stat(squeeze(u(ln,lt,d,:)),99);
%               [v_trend(ln,lt,d) v_plusminus(ln,lt,d), v_sig(ln,lt,d)] = trend_stat(squeeze(v(ln,lt,d,:)),99);
%               [w_trend(ln,lt,d) w_plusminus(ln,lt,d), w_sig(ln,lt,d)] = trend_stat(squeeze(w(ln,lt,d,:)),99);
% 
%          end
%      end
%  end
% 
% clear mon nmon d ln lt
%      
% % NaN out regions where trends are not significant
% u_trend(u_sig == 0) = NaN;
% v_trend(v_sig == 0) = NaN;
% w_trend(w_sig == 0) = NaN;
% 
% save EUC_Source_Waters.mat u_trend v_trend w_trend...
%     u_clim_trend v_clim_trend w_clim_trend...
%     u_clim v_clim w_clim u_anom v_anom w_anom...
%     lat lon depth 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Projection Time
%  
%  close all
%  
%  
% % set figure window size
% fig1 = figure('Color',[1 1 1]);
% ss=get(0,'ScreenSize');
% set(fig1,'Position',[ss(1)+15 ss(2)+25 ss(3) ss(4)-500])
% winsize = get(fig1,'Position');
% winsize(1:2) = [0 0];
% numframes=12;
% A=moviein(numframes,fig1,winsize);
% set(fig1,'NextPlot','replacechildren')
% month = {'January', 'February','March','April','May','June','July','August','September','October','November','December'}; 
% 

     
%% Indexing domain bounds
% NOTE: when implementing indices, use e.g. I_lon(1) and I_lon(end) for
% the data that make of the extremes edges of the domain of interest
 
 I_lon = find(lon>=160 & lon<=260);
 I_lat = find(lat<=4 & lat>=-4);     
 I_depth = find(depth>=0 & depth<=400);
 
 slon = lon(I_lon);
 slat = lat(I_lat);
 sdepth = depth(I_depth);

 nlon = length(slon);
 nlat = length(slat);
 
 
%%
% Trends (filled contours) w/ climatology (regular contours)
for month = 1:12

% CLIMATOLOGY; data are in dimensions: <lon x lat x depth x month>
% lon x lat planes: specify single depth
planeimg1(:,:,month) = squeeze(u_clim_trend(I_lon,I_lat,I_depth(1),month));        % upper plane
planeimg2(:,:,month) = squeeze(u_clim_trend(I_lon,I_lat,I_depth(end),month));      % lower plane

% lon x depth planes: specify single latitude
planeimg3(:,:,month) = squeeze(u_clim_trend(I_lon,I_lat(1),I_depth,month));          % southern plane
planeimg4(:,:,month) = squeeze(u_clim_trend(I_lon,I_lat(end),I_depth,month));        % northern plane

% lat x depth planes: specify single longitude
planeimg5(:,:,month) = squeeze(u_clim_trend(I_lon(1),I_lat,I_depth,month));          % western plane
planeimg6(:,:,month) = squeeze(u_clim_trend(I_lon(end),I_lat,I_depth,month));        % eastern plane


% TREND; data are in dimensions: <lon x lat x depth x month>
% lon x lat planes: specify single depth
planeimg7(:,:,month) = squeeze(u_clim_trend(I_lon,I_lat,I_depth(1),month));        % upper plane
planeimg8(:,:,month) = squeeze(u_clim_trend(I_lon,I_lat,I_depth(end),month));      % lower plane

% lon x depth planes: specify single latitude
planeimg9(:,:,month) = squeeze(u_clim_trend(I_lon,I_lat(1),I_depth,month));          % southern plane
planeimg10(:,:,month) = squeeze(u_clim_trend(I_lon,I_lat(end),I_depth,month));        % northern plane

% lat x depth planes: specify single longitude
planeimg11(:,:,month) = squeeze(u_clim_trend(I_lon(1),I_lat,I_depth,month));          % western plane
planeimg12(:,:,month) = squeeze(u_clim_trend(I_lon(end),I_lat,I_depth,month));        % eastern plane

end



%% Save To File that can be managaged by personal computer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% POJECTING DATA - LATEST VERSION ON PERSONAL COMPUTER

% % get the corners of the domain in which the data occurs.
% %Is this the best way to do this?
% 
% 
% min_x = min(min(slon));
% min_y = min(min(slat));
% max_x = max(max(slon));
% max_y = max(max(slat));
% 
% min_y3 = min(min(slat));
% max_y3 = max(max(slat));
% 
% min_y4= min(min(slat));
% max_y4 = max(max(slat));
% 
% min_y2 = max_y + 10;
% min_y2 = min_y2 +(max_y-min_y);
% 
% % desired z position of the image plane.
% imgzposition = 1100;          
% imgzposition1 = 500;
% imgzposition2 = 180;       
% 
% % latitudinal limits
% imgzposition3 = 0;          
% imgzposition4 = 4;
% 
% % longitudinal limits
% imgzposition5 = 140;
% imgzposition6 = 280;
% 
% %% lon x lat panels
% 
% % set inital background color to code NaN values in contour plots
% whitebg([0.1 0.1 0.1])  % dark grey background
% 
% % upper panel
%     surf([min_x max_x],[min_y max_y],repmat(imgzposition1, [2 2]),...
%     0*ones(size(planeimg3)),'facecolor','texture','FaceAlpha',1) 
%     
%     % trend contours
%     [C h] = contourf(slon,slat,planeimg1',100,'LineWidth',2);
%     set(h, 'EdgeColor','none')
%     
%     % climatolorgy contours
%     [C p] = contour(slon,slat,planeing1', 200, 'LineWidth',2);
% 
%     % panel rotation
%     zdir = [1 0 0];
%     origin1 = [0 0 imgzposition1];
%     rotate(get(h,'children'),zdir,90,origin1)
%     rotate(get(h,'children'),zdir,90,origin1)
% 
% hold on;
% 
% % lower panel
% surf([min_x max_x],[min_y max_y],repmat(-1*imgzposition+.5, [2 2]),...
%        0*ones(size(planeimg3)),'facecolor','texture','FaceAlpha',1)
%    
% %% lon x depth panels  
% 
% % Equator panel
% h12 = surf([min_x max_x],-1*[min(sdepth) max(sdepth)],repmat(-1*imgzposition3-.1, [2 2]),...
%        NaN(size(planeimg3)),'facecolor','texture','FaceAlpha',0);
%    
% zdir = [1 0 0];
% origin1 = [0 0 imgzposition3];
% rotate(h12,zdir,90,origin1)
% 
% % upper panel
% h12 = surf([min_x max_x],-1*[min(sdepth) max(sdepth)],repmat(-1*imgzposition4-.1, [2 2]),...
%        NaN(size(planeimg7)),'facecolor','texture','FaceAlpha',1);
%    
% [C h] = contourf(slon,slat,planeimg7',100,'LineWidth',2);
% set(h, 'EdgeColor','none')
%    
% zdir = [1 0 0];
% origin1 = [0 0 imgzposition3];
% rotate(h12,zdir,90,origin1)
% rotate(get(h,'children'),zdir,0,origin1) 
% 
% % lower panel
% h12 = surf([min_x max_x],-1*[min(sdepth) max(sdepth)],repmat(imgzposition4-.1, [2 2]),...
%        -.75*ones(size(planeimg3)),'facecolor','texture','FaceAlpha',1);
%    
% zdir = [1 0 0];
% origin1 = [0 0 imgzposition3];
% rotate(h12,zdir,90,origin1)
% 
% %% lat x depth panels   
% 
% % West Plane
% h12 = surf([min(sdepth) max(sdepth)],-1*[min_y max_y],repmat(imgzposition5-.1, [2 2]),...
%        .6*ones(size(planeimg3)),'facecolor','texture','FaceAlpha',1);
%    
% zdir = [0 1 0];
% origin1 = [0 0 imgzposition3];
% rotate(h12,zdir,90,origin1)
% 
% 
% % East Plane
% h12 = surf([min(sdepth) max(sdepth)],-1*[min_y max_y],repmat(imgzposition6-.1, [2 2]),...
%       .75*ones(size(planeimg3)),'facecolor','texture','FaceAlpha',1);
%    
% zdir = [0 1 0];
% origin1 = [0 0 imgzposition3];
% rotate(h12,zdir,90,origin1)
% 
% %% Rotate view
% view(13,25);
% whitebg([.8 .8 .8])
% 
% xlim([imgzposition5 imgzposition6])
% zlim([-1*imgzposition imgzposition1])
% colormap([0.1 0.1 0.1;jet(250)]);
% caxis([-1.01 1])
% 
% set(gcf,'Units','normal')
% set(gca,'Position',[.08 .08 0.87 0.91])
% 
% end