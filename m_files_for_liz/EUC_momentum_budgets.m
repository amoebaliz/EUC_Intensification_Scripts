%% Testing the Momentum budget of the EUC/determining relative importance of terms
%% Get SODA necessary data
close all
clear all

ncid = netcdf.open ('u_eqpac_SODA_extended.nc','NC_NOWRITE');               % zonal velocity
ncid_2 = netcdf.open ('v_eqpac_SODA_extended.nc','NC_NOWRITE');             % meridional velocity
ncid_3 = netcdf.open ('w_eqpac_SODA_extended.nc','NC_NOWRITE');             % vertical velocity
ncid_4 = netcdf.open ('taux_eqpac_SODA_extended_lat10.nc','NC_NOWRITE');    % winds
ncid_5 = netcdf.open ('ssh_eqpac_SODA_extended.nc','NC_NOWRITE');           % SSH
ncid_6 = netcdf.open ('temp_eqpac_SODA_extended.nc','NC_NOWRITE');  
ncid_7 = netcdf.open ('salt_eqpac_SODA_extended.nc','NC_NOWRITE');


varid_lon = netcdf.inqVarID(ncid,'lon');    % longitude
varid_lat = netcdf.inqVarID(ncid,'lat');    % latitude
varid_depth = netcdf.inqVarID(ncid,'depth');    % depth
varid_time = netcdf.inqVarID(ncid,'time');      % time

varid_u_euc = netcdf.inqVarID(ncid,'u');
varid_v_euc = netcdf.inqVarID(ncid_2,'v');
varid_w_euc = netcdf.inqVarID(ncid_3,'w');
varid_taux = netcdf.inqVarID(ncid_4,'taux');
varid_ssh = netcdf.inqVarID(ncid_5,'ssh');
varid_sst = netcdf.inqVarID(ncid_6,'temp');
varid_salt = netcdf.inqVarID(ncid_7,'salt');

    
lon = netcdf.getVar(ncid,varid_lon);
lat = netcdf.getVar(ncid,varid_lat);
depth = netcdf.getVar(ncid,varid_depth);
time = netcdf.getVar(ncid,varid_time);

u_EUC = netcdf.getVar(ncid,varid_u_euc);
v_EUC = netcdf.getVar(ncid_2,varid_v_euc);
w_EUC = netcdf.getVar(ncid_3,varid_w_euc);
taux = netcdf.getVar(ncid_4,varid_taux);
ssh = netcdf.getVar(ncid_5,varid_ssh);
temp = netcdf.getVar(ncid_6,varid_sst);
salt = netcdf.getVar(ncid_7,varid_salt);

u_EUC(u_EUC<-10000) = NaN;
v_EUC(v_EUC<-10000) = NaN;
w_EUC(w_EUC<-10000) = NaN;
taux(taux<-10000) = NaN;
ssh(ssh<-10000) = NaN;
temp(temp<-10000) = NaN;
salt(salt<-10000) = NaN;

%% initialize some values
close all
dens = zeros(size(temp));
for k = 1:size(temp,4)
for j = 1:size(temp,2)
    dens(:,j,:,k) = sw_dens(squeeze(squeeze(salt(:,j,:,k))),squeeze(squeeze(temp(:,j,:,k))),depth');
end
end
% [c h] = contourf(lon, depth, squeeze(squeeze(mean(mean(dens(:,10:11,:,:),2),4)))',5000);
% set(h,'EdgeColor','none')
% set(gca,'YDir','reverse')
%%
grav = 9.81;
delta_Lon = 111.320*1000; % coverted to meters
delta_Lat = 110.574*1000; % coverted to meters
delta_z = diff(depth);
for j = 1:(length(delta_z)-1)
    delta_z_2(j) = delta_z(j)+delta_z(j+1);
end

for j = 3:18
    delta_z_3(j-2) = depth(j+1)-depth(j-1);
end

mu = 1*10^-3;  % dynmic viscosity kg/ms
rho = 1028;   % kg/m^3
pressure_term = zeros(1,length(lon)-2);

%% Spatial indexing of latitude

% I_lat = find(lat<=0.5 & lat>=-0.5); 
 I_lat = find(lat<=2 & lat>=-2);
 slat = lat(I_lat);
 
 % upsilon/nu

%% Average values of all fields over time
avg_u = squeeze(mean(u_EUC,4));
avg_v = squeeze(mean(v_EUC,4));
avg_w = squeeze(mean(w_EUC,4));
avg_ssh = squeeze(mean(ssh,3));
avg_dens = squeeze(mean(dens,4));

%% Generate Climatology for all fields
%   want the climatology for all fields averaged on the equator 
%   (i.e. velocity and what not) BUT - what differs for zonal velocity - 
%   need climatolofy for 0.25 N/S, not just the average values.

u_EUC_2 = u_EUC(:,:,:,:);
u_EUC_3 = permute(u_EUC_2,[4 2 1 3]);

v_EUC_2 = v_EUC(:,:,:,:);
v_EUC_3 = permute(v_EUC_2,[4 2 1 3]);

w_EUC_2 = w_EUC(:,:,:,:);
w_EUC_3 = permute(w_EUC_2,[4 2 1 3]);

dens_2 = dens;
dens_3 = permute(dens_2,[4 2 1 3]);

u_EUC_clim = zeros(12, length(lat),length(lon),length(depth));
v_EUC_clim=u_EUC_clim;
w_EUC_clim=u_EUC_clim;
dens_clim = u_EUC_clim;

% climatology of current for both sides of the equator (for du/dy)
for d=1:length(depth)
    [u_EUC_clim(:,:,:,d), u_EUC_anom] = climanom(squeeze(u_EUC_3(:,:,:,d)));
    [v_EUC_clim(:,:,:,d), v_EUC_anom] = climanom(squeeze(v_EUC_3(:,:,:,d)));
    [w_EUC_clim(:,:,:,d), w_EUC_anom] = climanom(squeeze(w_EUC_3(:,:,:,d)));
    [dens_clim(:,:,:,d), dens_anom] = climanom(squeeze(dens_3(:,:,:,d)));
end


ssh_2 = permute(ssh,[3 2 1]);


[ssh_clim, ssh_anom] = climanom(ssh_2);

save('momentum.mat', 'avg_u','avg_v', 'avg_w', 'avg_ssh','avg_dens', 'u_EUC_clim',...
    'v_EUC_clim','w_EUC_clim','ssh_clim','lat','lon','depth','time', ...
    'grav','delta_Lon' , 'delta_Lat' ,'delta_z_2','delta_z_3', ...
    'mu', 'rho','I_lat')

% %%
% % on the equator (for u)
% % % u_EUC_4 = squeeze(mean(u_EUC(:,I_lat,:,:),2));
% % u_EUC_5 = permute(u_EUC_4,[3 2 1]);
% % [u_EUC_clim, u_EUC_anom] = climanom(u_EUC_5); 
%  
% % clim_check = (u_EUC_clim_1 + u_EUC_clim_2)/2;
% % check_ans = (clim_check == u_EUC_clim);
% 
% % meridional and vertical velocity (for v and w)
% % v_EUC_2 = squeeze(mean(v_EUC(:,I_lat,:,:),2));
% % v_EUC_3 = permute(v_EUC_2,[3 2 1]);
% % [v_EUC_clim, v_EUC_anom] = climanom(v_EUC_3);
% % 
% % w_EUC_2 = squeeze(mean(w_EUC(:,I_lat,:,:),2));
% % w_EUC_3 = permute(w_EUC_2,[3 2 1]);
% % [w_EUC_clim, w_EUC_anom] = climanom(w_EUC_3);
% % 
% % % winds and ssh climatologies 
% % taux_2 = permute(taux,[3 2 1]);
% % taux_3 = mean(taux_2(:,I_lat_2));
% % 
% % ssh_2 = permute(ssh,[3 2 1]);
% % ssh_3 = mean(ssh_2, 2);
% % 
% % [Taux_clim, Taux_anom] = climanom(taux_3);
% % [SSH_clim, SSH_anom] = climanom(ssh_3);
% 
% %% calculate all terms
% 
% for h = 1:12
%    for j = 2:(length(lon)-1)
%        for d = 2:(length(depth)-1)
% %% Calculate udu/dx
%         u_diff_x = u_EUC_clim(h,d,j+1)-u_EUC_clim(h,d,j-1);  
%         udu_term(j-1,d-1,h) = u_EUC_clim(h,d,j)*(u_diff_x/delta_Lon);
% 
% %% Calculate vdu/dy
%         u_diff_y = u_EUC_clim_1(h,d,j)-u_EUC_clim_2(h,d,j);
%         vdu_term(j-1,d-1,h) = v_EUC_clim(h,d,j)*(u_diff_y/delta_Lat);
% %% Calculate wdu/dz
%         u_diff_z = u_EUC_clim(h,d-1,j)-u_EUC_clim(h,d+1,j);
%         wdu_term(j-1,d-1,h) = w_EUC_clim(h,d,j)*(u_diff_z/delta_z_2(d-1));
% 
% %% Calculate Fx (windstress and propagated to depth)
%         % get du/dz for all depths before ending depth for look
%         dudz(j-1,d-1,h) = u_diff_z/delta_z_2(d-1);
%        end
%        for k = 2:17
%             Fx(j-1,k-1,h)=(1/rho)*mu*((dudz(j-1,k+1,h)-dudz(j-1,k-1,h))/delta_z_3(k-1));
%        end
% %% Calculate dP/dx
%         
% % the horizontal pressure gradient is identical everywhere within the
% % fluid; so it doesn't matter how deep you're looking
% 
% % this simplifies to: g*deltaZ/deltaX
% % g = gravitational constant (9.81 m s^-2)
% % deltaX (distance in longitude is 1 degree and assumed constant: 111.320 km)
% % deltaZ (difference in sea surface height) requires SSH at lon(i-1), (i+1)    
%         delta_SSH = SSH_clim(h,:,j-1)-SSH_clim(h,:,j+1); 
%         pressure_term(h,j-1) = grav*delta_SSH/delta_Lon;
%     end
% end
% %%
% close all
% cmin = -.0000003;
% cmax =  .0000003;
% 
% fig1 = figure;
% set(fig1,'Position',[100 100 1200 800]); 
% 
% for k = 1:12
%     
%     
% subplot(5,1,1)
% m = pcolor(repmat(lon(21:319)',[18,1]),repmat(depth(2:19),[1,299]),double(repmat(pressure_term(k,20:end),[18,1])));
% set(gca,'YDir','reverse')
% set(m, 'EdgeColor', 'none')
% caxis([cmin cmax])
% title(k)    
%     
% subplot(5,1,2)
% m = pcolor(repmat(lon(21:319)',[16,1]),repmat(depth(3:18),[1,299]),double(squeeze(Fx(20:end,:,k)))');
% set(gca,'YDir','reverse')
% set(m, 'EdgeColor', 'none')
% caxis([cmin*1000 cmax*1000])
%   
% 
% subplot(5,1,3)
% m = pcolor(repmat(lon(21:319)',[18,1]),repmat(depth(2:19),[1,299]),double(udu_term(20:end,:,k))');
% set(gca,'YDir','reverse')
% set(m, 'EdgeColor', 'none')
% caxis([cmin cmax])
% 
% subplot(5,1,4)
% p = pcolor(repmat(lon(21:319)',[18,1]),repmat(depth(2:19),[1,299]),double(vdu_term(20:end,:,k))');
% %contourf(lon(21:319),depth(2:19),wdu_term(20:end,:,k)')
% set(gca,'YDir','reverse')
% %cBar = colorbar;
% set(p, 'EdgeColor', 'none')
% caxis([cmin cmax])
% 
% subplot(5,1,5)
% t = pcolor(repmat(lon(21:319)',[18,1]),repmat(depth(2:19),[1,299]),double(wdu_term(20:end,:,k))');
% %contourf(lon(21:319),depth(2:19),wdu_term(20:end,:,k)')
% set(gca,'YDir','reverse')
% %cBar = colorbar;
% set(t, 'EdgeColor', 'none')
% caxis([cmin cmax])
% 
% 
% pause(2)
% end
% %%
% 
% close all
% fig2 = figure;
% L = pcolor(repmat(lon(2:319)',[12,1]),repmat((1:12)',[1,318]),pressure_term);
% set(L, 'EdgeColor', 'none')
% cmin = -0.0000006;
% cmax =  0.0000006;
% caxis([cmin cmax])