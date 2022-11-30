close all
clear all
clc

load mean_fields_E.mat
load redblue
I_lat = find(latx>=-0.5&latx<=0.5);
depth = 1:40;

subplot(2,1,1)
resid_plot = (dudt5);


[~,h]=contourf(lon,depth,10^7*(squeeze(mean(resid_plot(:,I_lat,:),2))'),100);
set(h,'EdgeColor','none')
hold on

% Y axis
set(gca,'YDir','reverse')
ylim([10 30])

set(gca,'TickDir','out','Ytick',[10 20 30],'fontsize', 14,'fontweight', 'bold')
% X axis

set(gca,'TickDir','out','Xtick',[])
xlim([150 270])

% C axis
cmax = 5;
hcb = colorbar; colormap(redblue); caxis([-1*cmax cmax]);

subplot(2,1,2)
resid_plot = (zpgf5-ududx5-vdudy5-wdudz5);


[~,h]=contourf(lon,depth,10^7*(squeeze(mean(resid_plot(:,I_lat,:),2))'),100);
set(h,'EdgeColor','none')
hold on

% Y axis
set(gca,'YDir','reverse')
ylim([10 30])
set(gca,'TickDir','out','Ytick',[10 20 30],'fontsize', 14,'fontweight', 'bold')
% X axis

set(gca,'TickDir','out','Xtick',[])
xlim([150 270])


% C axis
cmax = 5;
hcb = colorbar; colormap(redblue); caxis([-1*cmax cmax]);

%%

% Isopycnal zonal momentum budget with spatially varying vertical eddy
% viscosity for the equatorial Pacific Ocean in extended SODA

lat = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
lon = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
time = ncread('SODA_2.2.6_Trop_Pac_u.cdf','TIME');
depth = ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');
u = ncread('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN');


% eliminate missing value fill value (-9.99e+33)
u(u>1000) = NaN;

% Indexing spatial bounds
 I_lat = find(lat<=2 & lat>=-2);     
 I_lon = find(lon>=149 &lon<=271);
 I_dep = find(depth==depth);
 time0 = 61;
 
 slat = lat(I_lat);
 slon = lon(I_lon);
 sdep = depth(I_dep);

 nt = length(time);
 ndep = length(sdep);
 nlon = length(slon);
 nlat= length(slat);
 
% zonal velocities averaged over latitudes (and time)
u_euc_eq = squeeze(nanmean(u(I_lon,I_lat,I_dep,time0:end),2));

% statistical analyses
euc_trend = -500*ones(ndep,nlon);
euc_plusminus = euc_trend;
euc_sig = euc_trend;

%close all
for j = 1: nlon
    for k = 1:ndep
        
        [euc_trend(k,j),euc_plusminus(k,j),euc_sig(k,j)]=trend_stat(squeeze(u_euc_eq(j,k,:)),99);
   
    end
end

euc_trend(euc_sig<1) = NaN;
u = euc_trend;

%% define constants

dx=111320/2;
dy=110574/2;
dz=2;
rho0=1028;
% interpolate to uniform depth grid

depth(1)=5.0;
depth2=[5:dz:400]';
u2=zeros(length(depth2),size(u,2));

warning('off','all');
for x=1:size(u,2)
        u2(:,x)=interp1(depth,squeeze(u(:,x)),depth2);
end; clear x y
warning('on','all');
%%
u=u2; clear u2
% compute potential density

rhoth=NaN(size(temp));
for z=1:length(depth)
        rhoth(:,:,z)=sw_pden(squeeze(salt(:,:,z)),squeeze(temp(:,:,z)),depth(z),0);
end; clear z salt

for x=51:310
    for y=4:16
        z_tcline=min(depth(diff(temp(x,y,:)) == min(diff(temp(x,y,:)))));
        if( isfinite(z_tcline) && z_tcline ~= depth(1))
            z_deep=z_tcline+100;
            a=[depth(1) z_tcline z_deep];
            b=[Av_sfc Av_tcline Av_deep];
            c=depth(1:find(depth==z_deep));
            % d=spline(a,b,c);
            d=pchip(a,b,c);
            Av(x,y,1:find(depth==z_deep))=d;
            Av(x,y,depth>z_deep)=Av_deep;
        end
    end
end; clear x y z_tcline z_deep a b c d Av_sfc Av_tcline Av_deep 

% convert to isopycnal layers

isopyc_layer_bnds=1020:0.2:1028;

u_iso=NaN(length(lon),length(lat),length(isopyc_layer_bnds)-1);


for x=1:length(lon)
    for y=1:length(lat)
        for layer=1:length(isopyc_layer_bnds)-1
            zrange=depth(squeeze(rhoth(x,y,:))>=isopyc_layer_bnds(layer)&squeeze(rhoth(x,y,:))<=isopyc_layer_bnds(layer+1));
            if isfinite(zrange)
                u_iso(x,y,layer)=nanmean(u(x,y,depth>=zrange(1)&depth<=zrange(length(zrange))),3);
 
            end
        end
    end
end; clear x y layer zrange

clear rhoth rhoth_iso

u=u_iso; clear u_iso

