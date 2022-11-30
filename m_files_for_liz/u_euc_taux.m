close all
clear all

ncid = netcdf.open ('u_eqpac_SODA_extended.nc','NC_NOWRITE');               % zonal velocity
ncid_4 = netcdf.open ('taux_eqpac_SODA_extended_lat10.nc','NC_NOWRITE');    % winds

varid_lon = netcdf.inqVarID(ncid,'lon');    % longitude
varid_lat = netcdf.inqVarID(ncid,'lat');    % latitude
varid_depth = netcdf.inqVarID(ncid,'depth');    % depth
varid_time = netcdf.inqVarID(ncid,'time');      % time

varid_u_euc = netcdf.inqVarID(ncid,'u');
varid_taux = netcdf.inqVarID(ncid_4,'taux');

lon = netcdf.getVar(ncid,varid_lon);
lat = netcdf.getVar(ncid,varid_lat);
depth = netcdf.getVar(ncid,varid_depth);
time = netcdf.getVar(ncid,varid_time);

u_EUC = netcdf.getVar(ncid,varid_u_euc);
taux = netcdf.getVar(ncid_4,varid_taux);


u_EUC(u_EUC<-10000) = NaN;
taux(taux<-10000) = NaN;
%% Spatial indexing of latitude

 I_lat_taux = find(lat<=1 & lat>=-1); 
 I_lat_EUC = find(lat<=2 & lat>=-2);

 I_lon = find(lon<=270 & lon>=150); 
 
 I_lon_EUC = find(lon<=200 & lon>=160); 
 
 I_depth = find(depth<=500 & depth>=50); 
 nt = length(time);
 
 %%
 close all 

 u_max = zeros(1,nt);
 avg_taux = u_max;
 
 for n = 1:nt
     u_max(n) = max(max(max(u_EUC(I_lon_EUC,I_lat_EUC,I_depth,n))));
     avg_taux(n) = nansum(nansum(taux(I_lon,I_lat_taux,n)))/(numel(I_lat_taux)*numel(I_lon));
 end
 
 avg_taux_2 = runmean(avg_taux',13);
 u_max_2 = runmean(u_max',13);
 
 %plot(avg_taux_2, u_max_2,'bo')
 plot(avg_taux, u_max,'bo')
 
%[p, s] = polyfit(avg_taux_2(13:end-13), u_max_2(13:end-13), 2);
[p, s] = polyfit(avg_taux, u_max, 2);
Output = polyval(p,avg_taux);
Correlation = corrcoef(u_max, Output);


hold on
y = p(1)*avg_taux.^2 +p(2)*avg_taux+p(3);
plot(avg_taux,y,'k-')

begin_taux = mean(avg_taux(1:120));
begin_u_euc = mean(u_max(1:120));
taux_last = mean(avg_taux((end-120):end));
u_max_last = mean(u_max((end-120):end));

%begin_taux = mean(avg_taux_2(13:120));
%begin_u_euc = mean(u_max_2(13:120));
%taux_last = mean(avg_taux_2((end-140):end-13));
%u_max_last = mean(u_max_2((end-140):end-13));

plot(begin_taux,begin_u_euc,'r*')
plot(taux_last,u_max_last,'g*')