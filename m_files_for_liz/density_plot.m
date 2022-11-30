clear all
close all
clc

% SET DIFFERENCING SCHEME: 1 = central, 2 = consecutive

    diff_scheme = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lon=ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
lat=ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
depth=ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');
time=ncread('SODA_2.2.6_Trop_Pac_u.cdf','TIME');

temp=ncread('SODA_2.2.6_Trop_Pac_temp.cdf','TEMP_ENS_MN');
salt=ncread('SODA_2.2.6_Trop_Pac_salt.cdf','SALT_ENS_MN');

% Indexing spatial bounds
 I_lat = find(lat<=2 & lat>=-2);     
 I_lon = find(lon>=149 &lon<=271);
 I_dep = find(depth==depth);
 time0 = 61;

 slon = lon(I_lon);
 
[a b c d] = size(temp);


rho=NaN(length(I_lon),length(I_lat),length(I_dep),(d-time0+1));


for t = time0:length(time)
    for z=1:length(I_dep)
        rho(:,:,z,t-time0+1)=sw_dens(squeeze(salt(I_lon,I_lat,z,t)),squeeze(temp(I_lon,I_lat,z,t)),depth(z));
    end; clear z
end
clear salt temp

if diff_scheme == 1
    % vertical central differencing/ the distance (in meters)
    rho_dz = rho(:,:,1:end-2,:)-rho(:,:,3:end,:);
    dz = depth(3:end)-depth(1:end-2);
    dz_mat = repmat(permute(dz,[3 2 1]),[length(I_lon), length(I_lat),1,(d-time0+1)]);
    dens_grad = rho_dz./dz_mat;
    
    sdep = depth(2:end-1);

else
    % vertical differencing
    dz = depth(2:end)-depth(1:end-1);
    rho_dz = rho(:,:,1:end-1,:)-rho(:,:,2:end,:);
    dz_mat = repmat(permute(dz,[3 2 1]),[length(I_lon), length(I_lat),1,(d-time0+1)]);
    dens_grad = rho_dz./dz_mat;
    
    laywidth = zeros(length(depth)-1,1);
    
    for nlayer = 1:length(depth)
        if nlayer == 1
            laywidth(nlayer) = 2*depth(1);
        else
            laywidth(nlayer) = 2*(depth(nlayer)-sum(laywidth));
        end
    end
    
    sdep = depth(1:end-1)+(laywidth(1:end-1)/2);

end

%%%%

clear rho_dz dz_mat

%%%%

avg_grad = squeeze(mean(mean(dens_grad,4),2));
trend_grad = squeeze(mean(dens_grad,2));

dens_trend = NaN(length(slon),length(sdep));
dens_plusminus = dens_trend;
dens_sig = dens_trend;

for j = 1:length(slon)
    for k = 1:length(sdep)
        [dens_trend(j,k) dens_plusminus(j,k) dens_sig(j,k)] = trend_stat(squeeze(trend_grad(j,k,:)),99);
    end
end

%%%%%

save fig3_dens.mat slon sdep dens_trend dens_sig avg_grad

%% plotting data
close all
clc

figure('units','normal','outerposition',[.1 .5 .3 .6],'Color',[1 1 1])

load redblue
colormap(redblue)

dens_trend(dens_sig == 0) = NaN;

[a b] = contourf(slon,sdep,12*100*(10^3)*dens_trend',500);
set(b,'EdgeColor','none')
set(gca,'YDir','reverse')

cmin = -5;
cmax = 5;

caxis([cmin cmax])

hold on

[q,h] = contour(slon, sdep, avg_grad',[-.04,-.03,-.02,-.01],'k','LineWidth',2);
clabel(q,h, 'manual','fontsize', 13,'fontweight', 'bold','rotation',0)


%[q,h] = contour(slon,sdep,u_euc_avg',[-.40,-.20],'--b','LineWidth',2);
%[qq, hh] = contour(slon,sdep,u_euc_avg',[0, 0, 0],':k', 'LineWidth',2);
%[qqq,hhh] = contour(slon,sdep,u_euc_avg',[.20,.40,.60,.80],'k', 'LineWidth',2);

colorbar
