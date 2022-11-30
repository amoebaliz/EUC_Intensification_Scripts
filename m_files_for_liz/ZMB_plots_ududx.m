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

%% average ududx

load ududx.mat

avg_ududx = squeeze(mean(ududx(:,2,:,:),4));

clear ududx
%%
for j = 1:17
disp(j)
    % set time bounds for given subset 
    % Important because need array sizes to stay within MATLAB memory limits
   
   lon_bnds = 13*(j-1)+1:13*j;

load ududx.mat
ududx = ududx(lon_bnds,2,:,:);
A = squeeze(cumsum(ududx,4));

clear ududx

    if j == 1
    save A A
    else
        A_new = A; load A.mat; A = cat(1,A,A_new);
        save A A; clear A A_new lon_bnds

    end

end
%%
load A.mat

ududx_trend = NaN(size(A,2),size(A,1));
ududx_sig = ududx_trend;

for j = 1:size(A,1)
    for k = 1:size(A,2)
        [ududx_trend(k,j) ~, ududx_sig(k,j)] = trend(squeeze(A(j,k,:)),99);
    end
end

save ududx_trend ududx_trend ududx_sig slon depth2 avg_ududx

ududx_trend(ududx_sig == 0) = NaN;

clear A

slon = lon(I_lon);

%% plot lon x depth contour plots

close all

load redblue
colormap(redblue)
[a b] = contourf(slon(2:end-1),depth2',12*100*ududx_trend,500);
set(b,'EdgeColor','none')

set(gca, 'YDir', 'reverse')


cmax = 3*10^(-4);
cmin = -1*cmax;

caxis([cmin cmax])

colorbar


%%
hold on

[q, h] = contour(slon(2:end-1),depth2',avg_ududx',[0 0 0],'k','LineWidth',2);

[qq, hh] = contour(slon(2:end-1),depth2',avg_ududx',(.5:.5:1.5)*10^(-7),'k');
[qqq, hhh] = contour(slon(2:end-1),depth2',avg_ududx',(.5:.5:1.5)*-10^(-7),':k');



