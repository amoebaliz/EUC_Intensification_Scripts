close all; clear all; clc

%% define constants

global dx dy dz rho0 g Ah depth2

dx=111320/2;
dy=110574/2;
dz=2;
rho0=1028;
g=-9.81;
Ah=1.5*10^3; % Wallcroft et al. (2005, GRL)
% For Av, see constant vs. varying options below
depth2=(5:dz:400)';


%% load data from 2.2.6

global lon lat depth time

lon = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
lat = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
depth = ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');
time = ncread('SODA_2.2.6_Trop_Pac_u.cdf','TIME');

depth(1)=5.0;
%% Specify spatial/temporal bounds

global I_lon I_lat ntime nyearsub nyear ncycle

I_lon = find(lon >= 159.75 & lon <= 270.75);
I_lat = find(lat >= -0.75 & lat <= 0.75);

ntime = length(time);
nyearsub = 25;
nyear = ntime/12;
ncycle = ceil(nyear/nyearsub);

%% average zpgf

load vert_fric_dep.mat

avg_fricv = squeeze(mean(vert_fric_dep(:,:,:),3));

clear vert_fric_dep
%%
for j = 1:18
disp(j)
    % set time bounds for given subset 
    % Important because need array sizes to stay within MATLAB memory limits
   
   dep_bnds = 11*(j-1)+1:11*j;

load vert_fric_dep.mat
vert_fric_dep = vert_fric_dep(:,dep_bnds,:);
Fv = squeeze(cumsum(vert_fric_dep,3));

clear vert_fric_dep

    if j == 1
        save Fv Fv
    else
        Fv_new = Fv; load Fv.mat; Fv = cat(2,Fv,Fv_new);
        save Fv Fv; clear Fv Fv_new lon_bnds

    end

end
%%
load Fv.mat

fricv_trend = NaN(size(Fv,2),size(Fv,1));
fricv_sig = fricv_trend;

for j = 1:size(Fv,1)
    for k = 1:size(Fv,2)
        [fricv_trend(k,j) ~, fricv_sig(k,j)] = trend(squeeze(Fv(j,k,:)),99);
    end
end

%%
slon = lon(I_lon);
save fricv_trend fricv_trend fricv_sig slon depth2 avg_fricv

%fricv_trend(fricv_sig == 0) = NaN;

clear Fv

slon = lon(I_lon);

% plot lon x depth contour plots
%%
close all

load redblue
colormap(redblue)
[a b] = contourf(slon,depth2',12*100*fricv_trend,500);
set(b,'EdgeColor','none')

set(gca, 'YDir', 'reverse')


 cmax = 6 *10^(-4);
 cmin = -1*cmax;
% 
 caxis([cmin cmax])

colorbar


%%
hold on

[q, h] = contour(slon,depth2',avg_fricv',[0 0 0],'k','LineWidth',2);

[qq, hh] = contour(slon,depth2',avg_fricv',(.5:.5:2)*10^(-4),'k');
[qqq, hhh] = contour(slon,depth2',avg_fricv',(.5:.5:2)*-10^(-4),':k');



