close all
clear all

ncid = netcdf.open ('ssh_eqpac_SODA_extended.nc','NC_NOWRITE');           % SSH
ncid_2 = netcdf.open ('temp_eqpac_SODA_extended.nc','NC_NOWRITE');  
ncid_3 = netcdf.open ('salt_eqpac_SODA_extended.nc','NC_NOWRITE');

varid_lon = netcdf.inqVarID(ncid_2,'lon');        % longitude
varid_lat = netcdf.inqVarID(ncid_2,'lat');        % latitude
varid_depth = netcdf.inqVarID(ncid_2,'depth');    % depth

varid_ssh = netcdf.inqVarID(ncid,'ssh');
varid_sst = netcdf.inqVarID(ncid_2,'temp');
varid_salt = netcdf.inqVarID(ncid_3,'salt');

lon = netcdf.getVar(ncid_2,varid_lon);
lat = netcdf.getVar(ncid_2,varid_lat);
depth = netcdf.getVar(ncid_2,varid_depth);

ssh = netcdf.getVar(ncid,varid_ssh);
temp = netcdf.getVar(ncid_2,varid_sst);
salt = netcdf.getVar(ncid_3,varid_salt);

ssh(ssh<-10000) = NaN;
temp(temp<-10000) = NaN;
salt(salt<-10000) = NaN;

%% Constants
grav = 9.81;
delta_Lon = 111.320*1000; % coverted to meters
rho = 1028;   % kg/m^3
pressure_term = zeros(320,20);

%%
close all
dens = zeros(320,20);
for k = 800
    for j = 10
        dens = sw_dens(squeeze(squeeze(salt(:,j,:,k))),squeeze(squeeze(temp(:,j,:,k))),depth');
    end
end

subplot(3,1,1)
m = pcolor(repmat(lon',[20,1]),repmat(depth,[1,length(lon)]),double(dens'));
set(m, 'EdgeColor', 'none')
set(gca,'Ydir','reverse')
colorbar
title({'AUGUST 1937','Density Distribution (kg/m^3)'})

%% calculate depth intervals for pressure gradient term
y = zeros(length(depth)+1,1);
for j = 1:length(depth)-1
    y(j+1) = (depth(j+1)-depth(j))/2;
end
y(1) = depth(1);
y(end) = y(end-1);

depth_int = zeros(length(depth),1);
for j = 1:length(depth)
    depth_int(j) = y(j)+y(j+1);
end
%% Calculate Pressure Gradient Field

press_field_1 = zeros(320,20);

for j = 1:length(lon)
    %for k  = 1:length(lat)
        
        % density driven current:
        press_field_1(j,:) = dens(j,:).*depth_int';
        % adding on SSH component at the surface level:
        press_field_1(j,1) = press_field_1(j,1)+(dens(j,1)*ssh(j,10,800));
        
    %end
end



% Integrate pressure layers over depth
press_field_2 = grav*cumsum(press_field_1,2);

subplot(3,1,2)
m = pcolor(repmat(lon',[20,1]),repmat(depth,[1,length(lon)]),double(press_field_2'));
set(m, 'EdgeColor', 'none')
set(gca,'Ydir','reverse')
caxis([-5*10^6 5*10^6])
colorbar
title('Pressure Field (Pa)')

%% calculate all terms 

for j = 2:(length(lon)-1)
        
    %Calculate dP/dx
    pressure_term(j,:) = -1*(1/rho)*(press_field_2(j+1,:)-press_field_2(j-1,:))/delta_Lon;
    
end

subplot(3,1,3)
m = pcolor(repmat(lon',[20,1]),repmat(depth,[1,length(lon)]),pressure_term');
colorbar
set(m, 'EdgeColor', 'none')
set(gca,'Ydir','reverse')
caxis([-1*10^(-6) 1*10^-6])
title('Accleration due to Pressure Gradient (m/s^2)')
