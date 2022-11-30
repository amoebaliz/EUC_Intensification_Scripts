ncid_7 = netcdf.open ('prmsl.mon.mean.nc','NC_NOWRITE');

    % sea level pressure variable Ids. 
    varid_lat_prmsl = netcdf.inqVarID(ncid_7,'lat');
    varid_lon_prmsl = netcdf.inqVarID(ncid_7,'lon');
    varid_prmsl = netcdf.inqVarID(ncid_7,'prmsl');

 prmsl = netcdf.getVar(ncid_7,varid_prmsl);
 prmsl_lat = netcdf.getVar(ncid_7,varid_lat_prmsl);
 prmsl_lon = netcdf.getVar(ncid_7,varid_lon_prmsl);
 
 
 I_lat = find(prmsl_lat<=7 & prmsl_lat>=-7);      
 I_lon = find(prmsl_lon>=160 & prmsl_lon<=260);
 
 close all
for j = 1:12
    
    T = squeeze(prmsl(I_lon,I_lat,j)); 
    
    
    [C h] = contourf(prmsl_lon(I_lon),prmsl_lat(I_lat),T',100,'LineWidth',2);
    set(h, 'EdgeColor','none')
    caxis([-9200 -8000])
    colorbar
pause(.5)
end