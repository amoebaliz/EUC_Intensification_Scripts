

% This is a script to run the isopycnal zonal momentum budget on Isabela.
%
% Notes:
%
% You need to make sure the "seawater" library is in your path, because the
% ZMB code utilizes sw_dens and sw_pden (which themselves call other codes
% in that library). I acquired version 3.3 here from
% http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm
%
% The ZMB code requires input fields for just one time at a time. You can
% write a loop to loop through times and have it compute the ZMB terms for
% each time, saving the result of each iteration into a pre-defined array.
% 
% The example below works on Isabela to open a MAT file containing the SODA
% v2.2.6 record mean fields, pass them to the ZMB code, and plot the
% resulting 7 terms (u*du/dx, etc.) as well as the sum (du/dt).
%
% Because of all of the regridding steps and converting from z-levels to
% isopycnal (constant potential density) layers, it takes a couple of
% minutes for the code to run for just one time step.

close all, clear all
addpath('/data1/mcode/')
addpath('/data1/mcode/seawater_ver3_3/')
load redblue.mat

% load zonal velocity info to get dimension info

lon = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
lat = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
depth = ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');
time = ncread('SODA_2.2.6_Trop_Pac_u.cdf','TIME');

nmonths = length(time);
nyears = nmonths/12;

% preallocate space for saved variables

u_iso = NaN(320,19,40,nyears); 

dudt = u_iso; ududx = u_iso; vdudy = u_iso; wdudz = u_iso;
zpgf = u_iso; corf = u_iso; frich = u_iso; fricv = u_iso; depth_iso = u_iso;
%%

j2 = 0;
for j = 3:12:nmonths
    j2 = j2 +1;
    % load specific time field
    u = squeeze(ncread('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN', ...
        [1 1 1 j],[inf inf inf 1]));
    
    v = squeeze(ncread('SODA_2.2.6_Trop_Pac_v.cdf','V_ENS_MN', ...
        [1 1 1 j],[inf inf inf 1]));
    
    w = squeeze(ncread('SODA_2.2.6_Trop_Pac_w.cdf','W_ENS_MN', ...
        [1 1 1 j],[inf inf inf 1]));
    
    temp = squeeze(ncread('SODA_2.2.6_Trop_Pac_temp.cdf','TEMP_ENS_MN', ...
        [1 1 1 j],[inf inf inf 1]));
    
    salt = squeeze(ncread('SODA_2.2.6_Trop_Pac_salt.cdf', 'SALT_ENS_MN', ...
        [1 1 1 j],[inf inf inf 1]));
    
    ssh = squeeze(ncread('SODA_2.2.6_Trop_Pac_ssh.cdf','SSH_ENS_MN', ...
        [1 1 j],[inf inf 1]));
    
    
    [lon_iso, lat_iso, depth_iso(:,:,:,j2),u_iso(:,:,:,j2), dudt(:,:,:,j2), ududx(:,:,:,j2),...
        vdudy(:,:,:,j2), wdudz(:,:,:,j2), zpgf(:,:,:,j2), corf(:,:,:,j2),...
        frich(:,:,:,j2),fricv(:,:,:,j2)] = ...
        zmb_isopycnal(lon,lat,depth,u,v,w,temp,salt,ssh);
    
%     [lon, lat, depth,u_iso, dudt, ududx,...
%         vdudy, wdudz, zpgf, corf,...
%         frich,fricv] = ...
%         zmb_isopycnal(lon,lat,depth,u,v,w,temp,salt,ssh);

disp(j2)
end
%%


figure(1)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(ududx(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('u du/dx (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(2)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(vdudy(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('v du/dy (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(3)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(wdudz(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('w du/dz (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(4)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(zpgf(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('-\rho^{-1} dp/dx (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(5)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(corf(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('2 \Omega sin(\phi) (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(6)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(frich(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('A_h (d^2u/dx^2 + d^2u/dy^2) (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(7)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(fricv(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('d/dz(A_v*du/dz) (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(8)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(dudt(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('du/dt (sum) (m s^{-2})')
axis([140 280 1 size(depth,3)])
