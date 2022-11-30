close all 
clear all

% [x y z] = sphere(24);
% surf(x,y,z)

ncid_1 = netcdf.open ('u_eqpac_SODA_extended.nc','NC_NOWRITE');

% get variable IDs for zonal
    varid_time = netcdf.inqVarID(ncid_1,'time');    
    varid_lat = netcdf.inqVarID(ncid_1,'lat');
    varid_lon = netcdf.inqVarID(ncid_1,'lon');
    varid_depth = netcdf.inqVarID(ncid_1,'depth');
    varid_u = netcdf.inqVarID(ncid_1,'u');
    
%  get variables
    time = netcdf.getVar(ncid_1,varid_time);
    lat = netcdf.getVar(ncid_1,varid_lat);
    lon = netcdf.getVar(ncid_1,varid_lon);
    depth = netcdf.getVar(ncid_1,varid_depth);
    
    u_euc = netcdf.getVar(ncid_1,varid_u);
    u_euc(u_euc<-10000) = NaN;
 %%
 % Indexing bounds
 I_lat = find(lat<=3 & lat>=-3);      
 I_lon = find(lon>=120 &lon<=300);
 I_dep = find(depth>=30 & depth<=350);
 
 slat = lat(I_lat);
 slon = lon(I_lon);
 sdep = depth(I_dep);

 nt = length(time);
 nlon = length(slon);
 nlat= length(slat);
 
T = nanmean(nanmean(u_euc(I_lon,I_lat,I_dep,:),2),4);

    contour(slon, sdep, squeeze(T)')
    set(gca,'YDir','reverse')
    colorbar
    %%
euc_surf_velocity = 0.2;

T = squeeze(nanmean(u_euc(I_lon, I_lat, I_dep,:),4));

I = find(T>=euc_surf_velocity);
[x y z] = ind2sub([320 12 14],I);
%%
close all
plot3(slon(x),slat(y),sdep(z),'o')
set(gca,'ZDir','reverse')
%%

n = 1;
Iu = [];

%I2 = zeros(1,length(slon));
%I3 = zeros(1,length(slat));

for j = 1:length(slon)
    
    I2 = find (x ==j);
        if length(I2)>=1
        Iu(n:n+length(I2)-1) = I2;
        n=length(Iu)+1;
        end

    for k = 1:length(slat)
       % I3 = find (x ==j);
        

%        if x == j & y == k
%            
%            I2 = find(slon(x) == slon(j) & slat(y) == slat(k));

%            max_dep(n) = max(sdep(I2));       
%            min_dep(n) = min(sdep(I2));
       
%            n=n+1;
%        end
        
    end
end


%%

surf(SLON,SLAT,Z2);


