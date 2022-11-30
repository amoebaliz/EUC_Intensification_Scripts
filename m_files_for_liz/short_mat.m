close all 
clear all
clc

ncid = netcdf.open ('SODA_2.2.6_Trop_Pac_u.cdf','NC_NOWRITE');
ncid2 = netcdf.open ('SODA_2.2.6_Trop_Pac_v.cdf','NC_NOWRITE');
ncid3 = netcdf.open ('SODA_2.2.6_Trop_Pac_w.cdf','NC_NOWRITE');
ncid4 = netcdf.open ('SODA_2.2.6_Trop_Pac_salt.cdf','NC_NOWRITE');
ncid5 = netcdf.open ('SODA_2.2.6_Trop_Pac_temp.cdf','NC_NOWRITE');
ncid6 = netcdf.open ('SODA_2.2.6_Trop_Pac_ssh.cdf','NC_NOWRITE');
ncid7 = netcdf.open ('SODA_2.2.6_Trop_Pac_taux.cdf','NC_NOWRITE');

% get variable IDs
    varid_lat = netcdf.inqVarID(ncid,'LAT142_161');
    varid_lon = netcdf.inqVarID(ncid,'LON241_560');
    varid_time = netcdf.inqVarID(ncid,'TIME');
    varid_depth = netcdf.inqVarID(ncid,'DEPTH1_20');

    varid_u = netcdf.inqVarID(ncid,'U_ENS_MN');
    varid_v = netcdf.inqVarID(ncid2,'V_ENS_MN');
    varid_w = netcdf.inqVarID(ncid3,'W_ENS_MN');
    varid_salt = netcdf.inqVarID(ncid4,'SALT_ENS_MN');
    varid_temp = netcdf.inqVarID(ncid5,'TEMP_ENS_MN');
    varid_ssh = netcdf.inqVarID(ncid6,'SSH_ENS_MN');
    varid_taux = netcdf.inqVarID(ncid7,'TAUX_ENS_MN');
    
    
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
    taux = netcdf.getVar(ncid7,varid_taux);

% eliminate missing value fill value (-9.99e+33)
u(u<-10000) = NaN;
v(v<-10000) = NaN;
w(w<-10000) = NaN;
salt(salt<-10000) = NaN;
temp(temp<-10000) = NaN;
ssh(ssh<-10000) = NaN;
taux(taux<-10000) = NaN;

% Indexing spatial bounds

 I_lat = find(lat<=1 & lat>=-1);  
 I_lon = find(lon>=200 &lon<=210);
 I_depth = find(depth<= 200);
 
 
 nt = length(time);
 slon = lon(I_lon);
 
 u_short = squeeze(u(I_lon,I_lat,I_depth,600:620));
 v_short = squeeze(v(I_lon,I_lat,I_depth,600:620));
 w_short = squeeze(w(I_lon,I_lat,I_depth,600:620));
 salt_short = squeeze(salt(I_lon,I_lat,I_depth,600:620));
 temp_short = squeeze(temp(I_lon,I_lat,I_depth,600:620));
 
 ssh_short = squeeze(ssh(I_lon,I_lat,600:620));
 taux_short = squeeze(ssh(I_lon,I_lat,600:620));
 
 save ashort.mat u_short v_short w_short salt_short temp_short ssh_short taux_short depth slon
 