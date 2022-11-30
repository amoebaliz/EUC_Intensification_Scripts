% 
% 
% 
% load zpgf.mat
% 
% slon = lon(I_lon);
% %%
% close all
% for j=1:30
% contourf(slon(2:end-1),depth2',squeeze(zpgf(:,:,:,j+1000))',50)
% set(gca, 'YDir', 'reverse')
% 
% colorbar
% 
% caxis([-5 5]*10^(-7))
% pause(1)
% end

%%

close all; clc
load vert_fric.mat

load redblue
colormap(redblue)

avg_vert_fric = mean(fricv,3);
%%
slon = lon(I_lon);

for j = 1:200
     h = pcolor(repmat(slon',[36 1]),repmat((3:38)',[1 length(slon)]),squeeze(fricv(:,:,1200+j))');
     set(h, 'EdgeColor','none')
     set(gca, 'YDir', 'reverse')
     
     caxis([-7 7]*10^(-6))
     colorbar
     pause(0.5)
end
%%
 close all

figure(2)

h= pcolor(repmat(slon',[36 1]),repmat((3:38)',[1 length(slon)]),avg_vert_fric');
     set(h, 'EdgeColor','none')
     set(gca, 'YDir', 'reverse')
     colormap(redblue)
     caxis([-7 7]*10^(-6))
     colorbar
     title('mean vertical friction')

     figure(3)

h= pcolor(repmat(slon',[36 1]),repmat((3:38)',[1 length(slon)]),nanmean(fricv,3)');
     set(h, 'EdgeColor','none')
     set(gca, 'YDir', 'reverse')
     colormap(redblue)
     caxis([-7 7]*10^(-6))
     colorbar
    title('NaN-mean vertical friction')

 %%
 close all
 
 [v_clim v_anom] =climanom(permute(fricv,[3 1 2]));
 
 for k = 1:5
 for j = 1:12
 
 h= pcolor(repmat(slon',[36 1]),repmat((3:38)',[1 length(slon)]),squeeze(v_clim(j,:,:))');
 set(h, 'EdgeColor','none')
     set(gca, 'YDir', 'reverse')
     colormap(redblue)
     caxis([-7 7]*10^(-6))
     colorbar
     xlabel('Longitude')
     ylabel('Isopycnal')
 
 pause(0.5)
 end
 end
 
 
 
 
     
 
%%
% load redblue
% 
% u = ncread('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN');
% lon = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
% depth = ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');
% lat = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
% %%
% I_lat = find(lat <=2 & lat >=-2);
% I_lon = find(lon>=160 &lon<=260);
% slon = lon(I_lon);
% 
% u_first = squeeze(mean(mean(u(I_lon,I_lat,:,1:20),4),2));
% u_last = squeeze(mean(mean(u(I_lon,I_lat,:,end-19:end),4),2));
% 
% close all
% 
% cmax = 1; cmin = -1*cmax; 
% 
% subplot(3,1,1)
% 
% [c h] = contourf(slon, depth, u_first',500);
% set(h, 'EdgeColor','none')
% set(gca, 'YDir', 'reverse')
% 
% colormap(redblue)
% caxis([cmin cmax])
% colorbar
% 
% title('Average zonal velocity over first 20 yrs in SODA')
% ylabel('Depth (m)')
% 
% subplot(3,1,2)
% 
% [c h] = contourf(slon, depth, u_last',500);
% set(h, 'EdgeColor','none')
% set(gca, 'YDir', 'reverse')
% 
% colormap(redblue)
% caxis([cmin cmax])
% colorbar
% 
% title('Average zonal velocity over last 20 yrs in SODA')
% ylabel('Depth (m)')
% 
% subplot(3,1,3)
% 
% [c h] = contourf(slon, depth, (u_last - u_first)',500);
% set(h, 'EdgeColor','none')
% set(gca, 'YDir', 'reverse')
% 
% colormap(redblue)
% caxis([cmin cmax])
% colorbar
% 
% title('Difference between first and last 20 years')
% ylabel('Depth (m)')
% xlabel('Longitude')
% 
% 
