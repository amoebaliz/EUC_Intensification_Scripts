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

load zpgf.mat

avg_zpgf = squeeze(mean(zpgf(:,:,:),3));

clear zpgf
%%
for j = 1:18
disp(j)
    % set time bounds for given subset 
    % Important because need array sizes to stay within MATLAB memory limits
   
   dep_bnds = 11*(j-1)+1:11*j;

load zpgf.mat
zpgf = zpgf(:,dep_bnds,:);
E = squeeze(cumsum(zpgf,3));

clear zpgf

    if j == 1
        save E E
    else
        E_new = E; load E.mat; E = cat(2,E,E_new);
        save E E; clear E E_new lon_bnds

    end

end
%%
load E.mat

zpgf_trend = NaN(size(E,2),size(E,1));
zpgf_sig = zpgf_trend;

for j = 1:size(E,1)
    for k = 1:size(E,2)
        [zpgf_trend(k,j) ~, zpgf_sig(k,j)] = trend(squeeze(E(j,k,:)),99);
    end
end

%%
slon = lon(I_lon);
save zpgf_trend zpgf_trend zpgf_sig slon depth2 avg_zpgf

zpgf_trend(zpgf_sig == 0) = NaN;

clear E

slon = lon(I_lon);

% plot lon x depth contour plots

close all

load redblue
colormap(redblue)
[a b] = contourf(slon(2:end-1),depth2',12*100*zpgf_trend,500);
set(b,'EdgeColor','none')

set(gca, 'YDir', 'reverse')


 cmax = 6 *10^(-4);
 cmin = -1*cmax;
% 
 caxis([cmin cmax])

colorbar


%%
hold on

[q, h] = contour(slon(2:end-1),depth2',avg_zpgf',[0 0 0],'k','LineWidth',2);

[qq, hh] = contour(slon(2:end-1),depth2',avg_zpgf',(.5:.5:2)*10^(-7),'k');
[qqq, hhh] = contour(slon(2:end-1),depth2',avg_zpgf',(.5:.5:2)*-10^(-7),':k');



