close all 
clear all

% [x y z] = sphere(24);
% surf(x,y,z)

ncid = netcdf.open ('third_SODA_subset.cdf','NC_NOWRITE');
varid_u = netcdf.inqVarID(ncid,'U_ENS_MN');
u = netcdf.getVar(ncid,varid_u);
u(u<-10000) = NaN;

ncid_1 = netcdf.open ('eqpac_SODA_2.2.6_ens.cdf','NC_NOWRITE');

% get variable IDs for zonal
   
    varid_lat = netcdf.inqVarID(ncid_1,'LAT1_300');
    varid_lon = netcdf.inqVarID(ncid_1,'LON');
    varid_depth = netcdf.inqVarID(ncid_1,'DEPTH');
    varid_u = netcdf.inqVarID(ncid_1,'U_ENS_MN');
    
%  get variables
 
    lat = netcdf.getVar(ncid_1,varid_lat);
    lon = netcdf.getVar(ncid_1,varid_lon);
    depth = netcdf.getVar(ncid_1,varid_depth);
    
    u_euc = netcdf.getVar(ncid_1,varid_u);
    u_euc(u_euc<-10000) = NaN;
    
    subplot(3,1,1)
    m = pcolor(repmat(lon',[length(lat),1]),repmat(lat,[1,length(lon)]),double(squeeze(u_euc(1,:,:,1))));
    set(m, 'EdgeColor', 'none')

    subplot(3,1,2)
    m = pcolor(repmat(lon',[length(lat),1]),repmat(lat,[1,length(lon)]),double(squeeze(u_euc(1,:,:,2))));
    set(m, 'EdgeColor', 'none')
    
    subplot(3,1,3)
    m = pcolor(repmat(lon',[length(lat),1]),repmat(lat,[1,length(lon)]),double(squeeze(u(1,:,:,:))));
    set(m, 'EdgeColor', 'none')
    
    