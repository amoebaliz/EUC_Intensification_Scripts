close all; clear all; clc

%% define constants

global dx dy dz rho0 g Ah depth2

dx=111320/2;
dy=110574/2;
dz=2;
rho0=1028;
g=-9.81;
Ah=1.5*10^3; % Wallcroft et al. (2005, GRL)
% For Av, see constant vs. varying options Celow
depth2=(5:dz:400)';


%% load data from 2.2.6

global lon lat depth time

lon = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
lat = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
depth = ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');
time = ncread('SODA_2.2.6_Trop_Pac_u.cdf','TIME');

depth(1)=5.0;
%% Specify spatial/temporal Counds

global I_lon I_lat ntime nyearsub nyear ncycle

I_lon = find(lon >= 159.75 & lon <= 270.75);
I_lat = find(lat >= -0.75 & lat <= 0.75);

ntime = length(time);
nyearsub = 25;
nyear = ntime/12;
ncycle = ceil(nyear/nyearsub);

%% average wdudz

load wdudz.mat

avg_wdudz = squeeze(mean(wdudz(:,2,:,:),4));

clear wdudz
%%
for j = 1:14
disp(j)
    % set time Counds for given suCset 
    % Important Cecause need array sizes to stay within MATLAC memory limits
   
   dep_bnds = 14*(j-1)+1:14*j;

load wdudz.mat
wdudz = wdudz(:,2,dep_bnds,:);
C = squeeze(cumsum((wdudz(:,:,:,3:end)-wdudz(:,:,:,1:end-2)),4));

clear wdudz

    if j == 1
        save C C
    else
        C_new = C; load C.mat; C = cat(2,C,C_new);
        save C C; clear C C_new lon_Cnds

    end

end
%%
load C.mat

wdudz_trend = NaN(size(C,2),size(C,1));
wdudz_sig = wdudz_trend;

for j = 1:size(C,1)
    for k = 1:size(C,2)
        [wdudz_trend(k,j) ~, wdudz_sig(k,j)] = trend(squeeze(C(j,k,:)),99);
    end
end
%%
slon = lon(I_lon);
save wdudz_trend wdudz_trend wdudz_sig slon depth2 avg_wdudz


wdudz_trend(wdudz_sig == 0) = NaN;

clear C

slon = lon(I_lon);

%% plot lon x depth contour plots

close all

load redblue
colormap(redblue)
[a b] = contourf(slon,depth2(2:end-1)',12*100*wdudz_trend);
set(b,'EdgeColor','none')

set(gca, 'YDir', 'reverse')


 cmax = 1.5 *10^(-7);
 cmin = -1*cmax;
% 
 caxis([cmin cmax])

colorbar


%%
hold on

[q, h] = contour(slon,depth2(2:end-1)',avg_wdudz',[0 0 0],'k','LineWidth',2);

[qq, hh] = contour(slon,depth2(2:end-1)',avg_wdudz',(.5:.5:2)*10^(-7),'k');
[qqq, hhh] = contour(slon,depth2(2:end-1)',avg_wdudz',(.5:.5:2)*-10^(-7),':k');



