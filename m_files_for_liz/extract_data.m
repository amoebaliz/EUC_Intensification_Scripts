% SOvDA Scripts for SEC transports and taux
 
 close all
 clear all
 
    ncid = netcdf.open ('u_eqpac_SODA_extended_lat10.nc','NC_NOWRITE');
    ncid_2 = netcdf.open ('taux_eqpac_SODA_extended_lat10.nc','NC_NOWRITE');    
    % for looking at EUC 
    ncid_3 = netcdf.open ('u_eqpac_SODA_extended.nc','NC_NOWRITE'); 
    
    % Temperature and sea surface height data
    ncid_4 = netcdf.open ('sst_eqpac_SODA_extended.nc','NC_NOWRITE');
    ncid_5 = netcdf.open ('ssh_eqpac_SODA_extended.nc','NC_NOWRITE');
    ncid_6 = netcdf.open ('temp_eqpac_SODA_extended.nc','NC_NOWRITE');
    
    %global mean pressure 
    ncid_7 = netcdf.open ('prmsl.mon.mean.nc','NC_NOWRITE');
    
    % get variable IDs for SEC
    varid_lat_1 = netcdf.inqVarID(ncid,'lat');
    varid_lon = netcdf.inqVarID(ncid,'lon');
    varid_time = netcdf.inqVarID(ncid,'time');
    varid_u_sec = netcdf.inqVarID(ncid,'u');
    
    % get variable IDs for taux
    varid_taux = netcdf.inqVarID(ncid_2,'taux');
    
    % get variable IDs for EUC
    varid_u_euc = netcdf.inqVarID(ncid,'u');
    varid_depth_euc = netcdf.inqVarID(ncid,'depth');
    varid_lat_2 = netcdf.inqVarID(ncid,'lat');
    
    % temp variable IDs
    varid_sst = netcdf.inqVarID(ncid_4,'temp');
    varid_ssh = netcdf.inqVarID(ncid_5,'ssh');
    varid_ssh_lat = netcdf.inqVarID(ncid_5,'lat');
    varid_temp = netcdf.inqVarID(ncid_6,'temp');
    
    % sea level pressure variable Ids. 
    varid_lat_prmsl = netcdf.inqVarID(ncid_7,'lat');
    varid_lon_prmsl = netcdf.inqVarID(ncid_7,'lon');
    varid_prmsl = netcdf.inqVarID(ncid_7,'prmsl');

    % get variables (general)
    time = netcdf.getVar(ncid,varid_time);
    sec_lat = netcdf.getVar(ncid,varid_lat_1);
    euc_lat = netcdf.getVar(ncid_3,varid_lat_2);
    ssh_lat = netcdf.getVar(ncid_5,varid_ssh_lat);
    lon = netcdf.getVar(ncid,varid_lon);
    
    % variables for SEC (lon x lat x depth x time)
    u_sec = netcdf.getVar(ncid,varid_u_sec);
    
    % variables for taux (lon x lat x time)
    taux = netcdf.getVar(ncid_2,varid_taux);
    
    % variables for EUC (lon x lat x depth x time)
    u_euc = netcdf.getVar(ncid,varid_u_euc);
    depth = netcdf.getVar(ncid,varid_depth_euc);
    
    % variables for SST
    sst = netcdf.getVar(ncid_4,varid_sst);
    ssh = netcdf.getVar(ncid_5,varid_ssh);
    ssh_lat = netcdf.getVar(ncid_5,varid_ssh_lat);
    temp = netcdf.getVar(ncid_6,varid_temp);
    
    % variables for pressure and whatnot
    prmsl = netcdf.getVar(ncid_7,varid_prmsl);
    prmsl_lat = netcdf.getVar(ncid_7,varid_lat_prmsl);
    prmsl_lon = netcdf.getVar(ncid_7,varid_lon_prmsl);
 
 % eliminate missing value fill value (-9.99e+33)
 u_sec(u_sec<-10000) = NaN;
 taux(taux<-10000) = NaN;
 u_euc(u_euc<-10000) = NaN;
 

 %% Indexing spatial bounds
 
 depth_sec_2 = 15.0700;
 
 I_lat_1 = find(sec_lat<=7 & sec_lat>=-7);      % For SEC, Taux, SST, SSH
 I_lat_2 = find(euc_lat<=3 & euc_lat>=-3);      % For EUC, temp
 I_lat_3 = find(prmsl_lat<=7 & prmsl_lat>=-7);  % For SLP
 
 I_lon = find(lon>=160 &lon<=260);
 I_lon_2 = find(prmsl_lon>=160 & prmsl_lon<=260); 
 
 slat_1 = sec_lat(I_lat_1);
 slat_2 = euc_lat(I_lat_2);
 slat_3 = prmsl_lat(I_lat_3);
 
 slon = lon(I_lon);
 SLP_lon = prmsl_lon(I_lon_2);
 nt = length(time);
 nlon = length(slon);
 nlat_1 = length(slat_1);
 nlat_2 = length(slat_2);
 
%% Putting time in datenum form

time_2 = zeros(1,length(time));
 
 timen = time-.5;
 for j = 1:length(time)
    time_2(j) = addtodate(datenum(1960,1,1),double(timen(j)),'month');
 end
 
 %%
 
U_EUC = squeeze(mean(u_euc(I_lon,I_lat_2,:,:),2));
TEMP = squeeze(mean(temp(I_lon,I_lat_2,:,:),2));

% lat lon grids
U_SEC = squeeze(u_sec(I_lon,I_lat_1,:,:)); 
TAUX = squeeze(taux(I_lon,I_lat_1,:));
SSH = squeeze(ssh(I_lon,:,:));
SST = squeeze(sst(I_lon,:,:));

% SLP stuff
SLP = squeeze(prmsl(I_lon_2,I_lat_3,:));



save('EUC_data.mat','U_EUC', 'TEMP','U_SEC', 'TAUX','SSH', 'SST', 'SLP','slon', 'SLP_lon','slat_1', 'slat_3', 'depth', 'ssh_lat')





