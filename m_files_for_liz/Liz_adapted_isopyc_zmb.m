% This is a script to run the isopycnal zonal momentum budget on Isabela.
%
% Notes:
%
% You need to make sure the "seawater" library is in your path, because the
% ZMB code utilizes sw_dens and sw_pden (which themselves call other codes
% in that library). I acquired version 3.3 here from
% http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm
%
% The ZMB code requires input fields for just one time at a time. You can
% write a loop to loop through times and have it compute the ZMB terms for
% each time, saving the result of each iteration into a pre-defined array.
% 
% The example below works on Isabela to open a MAT file containing the SODA
% v2.2.6 record mean fields, pass them to the ZMB code, and plot the
% resulting 7 terms (u*du/dx, etc.) as well as the sum (du/dt).
%
% Because of all of the regridding steps and converting from z-levels to
% isopycnal (constant potential density) layers, it takes a couple of
% minutes for the code to run for just one time step.

close all; clear all; clc
addpath('/data1/mcode/')
addpath('/data1/mcode/seawater_ver3_3/')
load redblue.mat

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
nyearsub = 13;
nyear = ntime/12;
ncycle = ceil(nyear/nyearsub);

%% load and interpolate variables

du_cor_terms('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN', ...
   'SODA_2.2.6_Trop_Pac_v.cdf','V_ENS_MN',...
   'SODA_2.2.6_Trop_Pac_w.cdf','W_ENS_MN')

%%
pres_term('SODA_2.2.6_Trop_Pac_salt.cdf', 'SALT_ENS_MN',...
    'SODA_2.2.6_Trop_Pac_temp.cdf','TEMP_ENS_MN',...
    'SODA_2.2.6_Trop_Pac_ssh.cdf','SSH_ENS_MN')

%%
fric_terms('SODA_2.2.6_Trop_Pac_temp.cdf','TEMP_ENS_MN',...
    'SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN')

%% Friction terms

% Constant Av options

% Av=ones(size(u))*1.66*10^(-3); % Bryden & Brady (1985)
% Av=ones(size(u))*1*10^(-4); % MITgcm manual

% VERTICALLY VARYING COEFFICIENT OF VERTICAL EDDY VISCOSITY (Av) SEE QIAO
% AND WEISBERG (1997, JPO, FIG. 14) THERE ARE TWO OPTIONS: TIE INFLECTION
% POINT TO THE THERMOCLINE (MAX DT/DZ) OR TO THE EUC (MAX U)

Av_sfc=45*10^(-4); Av_tcline=3*10^(-4); Av_deep=15*10^(-4);

% OPTION 1: TIE INFLECTION POINT TO THE THERMOCLINE

% deep is depth of thermocline + 100 m. note that QW97 says the values
% below the EUC might be at least an order of magnitude smaller than this
% based on measurements.

temp = ncread('SODA_2.2.6_Trop_Pac_temp.cdf','TEMP_ENS_MN',...
  [I_lon(1) I_lat(1) 1 1],[length(I_lon) length(I_lat) inf inf]);
 temp(temp<-10000) = NaN;
% 
temp2 = interp1(depth, permute(temp,[3 2 1 4]), depth2);
temp3 = permute([mean(temp2(:,1:2,:,:),2) mean(temp2(:,2:3,:,:),2) mean(temp2(:,3:4,:,:),2)],[3 2 1 4]);
temp=temp3; clear temp2 temp3
%%
z_tcline = zeros(size(temp));

min(depth2(diff(temp) == min(diff(temp))));






%%










[lon,lat,depth,u,dudt,frich,fricv] = liz_zmb_isopycnal(lon,lat,depth,u,v,w,temp,salt,ssh);


% choose to work on the equator
I_lat = find(lat == 0);

% integrate dudt over time
dudt_integral = cumsum(dudt,4);

% refress dudt for all lat and all lon on equator

zmb_u_trend = NaN(length(lon)-2,length(depth)-2);
zmb_plusminus = zmb_u_trend;
zmb_sig = zmb_u_trend;

for j = 2: length(lon)-1
    for k = 3: length(depth)-2
        [zmb_u_trend(j,k), zmb_plusminus(j,k), zmb_sig(j,k)] = trend_stat(squeeze(dudt_integral(j,I_lat,k,:)),99);
    end
end

avg_u = squeeze(mean(mean(u(2:end-1,lat>=-0.5&lat<=0.5,3:end-2,:),2),4));

u_trend = NaN(length(lon), length(depth));
plusminus = u_trend;
sig = u_trend; 

for j = 1: length(lon)
    for k = 1: length(depth)
        [u_trend(j,k), plusminus(j,k), sig(j,k)] = trend_stat(squeeze(u(j,lat>=-0.5&lat<=0.5,k,:)),99);
    end
end

sublot(2,1,1)

[~,h]=contourf(lon(2:end-1),depth(3:end-2),zmb_u_trend',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); % caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; 
q = contour(lon(2:end-1),depth(3:end-2),avg_u',[-40,-20],'-.b','LineWidth',2);
qq = contour(lon(2:end-1),depth(3:end-2),[0, 0, 0],':k', 'LineWidth',2);
qqq = contour(lon(2:end-1),depth(3:end-2),avg_u',[20,40,60,80],'k', 'LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Depth'); title('regression of momentum balance terms (m s^{-2})')

subplot(2,1,2)

[~,h]=contourf(lon(2:end-1),depth(3:end-2),u_trend',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); % caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; 
q = contour(lon(2:end-1),depth(3:end-2),avg_u',[-40,-20],'-.b','LineWidth',2);
qq = contour(lon(2:end-1),depth(3:end-2),avg_u',[0, 0, 0],':k', 'LineWidth',2);
qqq = contour(lon(2:end-1),depth(3:end-2),avg_u',[20,40,60,80],'k', 'LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Depth'); title('regression of momentum balance terms (m s^{-2})')