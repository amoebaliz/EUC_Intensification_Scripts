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

%load record_mean_fields_v226.mat

%% load data from 2.2.6

ncid = netcdf.open ('SODA_2.2.6_Trop_Pac_u.cdf','NC_NOWRITE');
ncid2 = netcdf.open ('SODA_2.2.6_Trop_Pac_v.cdf','NC_NOWRITE');
ncid3 = netcdf.open ('SODA_2.2.6_Trop_Pac_w.cdf','NC_NOWRITE');
ncid4 = netcdf.open ('SODA_2.2.6_Trop_Pac_salt.cdf','NC_NOWRITE');
ncid5 = netcdf.open ('SODA_2.2.6_Trop_Pac_temp.cdf','NC_NOWRITE');
ncid6 = netcdf.open ('SODA_2.2.6_Trop_Pac_ssh.cdf','NC_NOWRITE');


% get variable IDs
    varid_lat = netcdf.inqVarID(ncid,'LAT142_161');
    varid_lon = netcdf.inqVarID(ncid,'LON241_560');
    varid_depth = netcdf.inqVarID(ncid,'DEPTH1_20');
    varid_time = netcdf.inqVarID(ncid,'TIME');

    varid_u = netcdf.inqVarID(ncid,'U_ENS_MN');
    varid_v = netcdf.inqVarID(ncid2,'V_ENS_MN');
    varid_w = netcdf.inqVarID(ncid3,'W_ENS_MN');
    varid_salt = netcdf.inqVarID(ncid4,'SALT_ENS_MN');
    varid_temp = netcdf.inqVarID(ncid5,'TEMP_ENS_MN');
    varid_ssh = netcdf.inqVarID(ncid6,'SSH_ENS_MN');

    
    
% get variables (general)
    time = netcdf.getVar(ncid,varid_time);
    lat = netcdf.getVar(ncid,varid_lat);
    lon = netcdf.getVar(ncid,varid_lon);
    depth = netcdf.getVar(ncid,varid_depth);
    
    u = netcdf.getVar(ncid,varid_u);
    v = netcdf.getVar(ncid2,varid_v);
    w = netcdf.getVar(ncid3,varid_w);
    salt = netcdf.getVar(ncid4,varid_salt);
    temp = netcdf.getVar(ncid5,varid_temp);
    ssh = netcdf.getVar(ncid6,varid_ssh);

% eliminate missing value fill value (-9.99e+33)
u(u<-10000) = NaN;
v(v<-10000) = NaN;
w(w<-10000) = NaN;
salt(salt<-10000) = NaN;
temp(temp<-10000) = NaN;
ssh(ssh<-10000) = NaN;


%%





[lon,lat,depth,u,dudt,ududx,vdudy,wdudz,zpgf,corf,frich,fricv] = zmb_isopycnal(lon,lat,depth,u,v,w,temp,salt,ssh);

figure(1)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(ududx(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('u du/dx (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(2)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(vdudy(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('v du/dy (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(3)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(wdudz(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('w du/dz (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(4)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(zpgf(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('-\rho^{-1} dp/dx (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(5)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(corf(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('2 \Omega sin(\phi) (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(6)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(frich(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('A_h (d^2u/dx^2 + d^2u/dy^2) (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(7)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(fricv(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('d/dz(A_v*du/dz) (m s^{-2})')
axis([140 280 1 size(depth,3)])

figure(8)
[~,h]=contourf(lon,1:size(depth,3),squeeze(mean(dudt(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('du/dt (sum) (m s^{-2})')
axis([140 280 1 size(depth,3)])
