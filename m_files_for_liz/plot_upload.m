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

%%

% save june_isopyc ...
%     ...
%     lon_iso lat_iso depth_iso u_iso dudt ududx vdudy wdudz zpgf corf

%%
% plotting the trends and the average u_iso for each of the saved fields

I_lat = find(lat<=0.5 & lat>=-.5);

% allocate space in advance
u_iso_trend = NaN(40,length(lon_iso));

u_iso_plusminus = u_iso_trend; u_iso_sig = u_iso_trend; 
ududx_trend = u_iso_trend; ududx_plusminus = u_iso_trend; ududx_sig = u_iso_trend;
vdudy_trend = u_iso_trend; vdudy_plusminus = u_iso_trend; vdudy_sig = u_iso_trend;
wdudz_trend = u_iso_trend; wdudz_plusminus = u_iso_trend; wdudz_sig = u_iso_trend;

zpgf_trend = u_iso_trend; zpgf_plusminus = u_iso_trend; zpgf_sig = u_iso_trend;
corf_trend = u_iso_trend; corf_plusminus = u_iso_trend; corf_sig = u_iso_trend;

fricv_trend = u_iso_trend; fricv_plusminus = u_iso_trend; fricv_sig = u_iso_trend;
frich_trend = u_iso_trend; frich_plusminus = u_iso_trend; frich_sig = u_iso_trend;

for j = 1:length(lon_iso)
    for k = 1:40
        
       [u_iso_trend(k,j),u_iso_plusminus(k,j),u_iso_sig(k,j)] = trend(squeeze(mean(u_iso(j,I_lat,k,:),2)),99);
       
       [ududx_trend(k,j),ududx_plusminus(k,j),ududx_sig(k,j)] = trend(squeeze(mean(ududx(j,I_lat,k,:),2)),99);
       [vdudy_trend(k,j),vdudy_plusminus(k,j),vdudy_sig(k,j)] = trend(squeeze(mean(vdudy(j,I_lat,k,:),2)),99);
       [wdudz_trend(k,j),wdudz_plusminus(k,j),wdudz_sig(k,j)] = trend(squeeze(mean(wdudz(j,I_lat,k,:),2)),99);
       
       [zpgf_trend(k,j),zpgf_plusminus(k,j),zpgf_sig(k,j)] = trend(squeeze(mean(zpgf(j,I_lat,k,:),2)),99);
       [corf_trend(k,j),corf_plusminus(k,j),corf_sig(k,j)] = trend(squeeze(mean(corf(j,I_lat,k,:),2)),99);
       
       [fricv_trend(k,j),fricv_plusminus(k,j),fricv_sig(k,j)] = trend(squeeze(mean(fricv(j,I_lat,k,:),2)),99);
       [frich_trend(k,j),frich_plusminus(k,j),frich_sig(k,j)] = trend(squeeze(mean(frich(j,I_lat,k,:),2)),99);
        
    end
end
%%
figure(2)
close all
subplot(2,1,2)

u_iso_trend2 = u_iso_trend;
u_iso_trend2(u_iso_sig==0)=NaN;

[~,h]=contourf(lon,1:size(depth_iso,3),u_iso_trend2,100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-3 5e-3]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('trend in u (m s^{-1}y^{-1})')
axis([140 280 1 size(depth_iso,3)])

delta_iso_u = squeeze(nanmean(nanmean(u_iso(:,I_lat,:,end-29:end),4),2))- squeeze(nanmean(nanmean(u_iso(:,I_lat,:,1:30),4),2));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_u',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-1 5e-1]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta u (m s^{-1})')
axis([140 280 1 size(depth_iso,3)])

%% ududx
close all
figure(2)
subplot(2,1,2)

ududx_trend2 = ududx_trend;
ududx_trend2(ududx_sig==0)=NaN;

[~,h]=contourf(lon,1:size(depth_iso,3),ududx_trend2,100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-1e-8 1e-8]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('trend in u du/dx (m s^{-2}y^{-1})')
axis([140 280 1 size(depth_iso,3)])


delta_ududx = squeeze(nanmean(nanmean(ududx(:,I_lat,:,end-29:end),4),2))- squeeze(nanmean(nanmean(ududx(:,I_lat,:,1:30),4),2));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_ududx',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta u du/dx (m s^{-2})')
axis([140 280 1 size(depth_iso,3)])

%%
 close all

figure(3)

vdudy_trend2 = vdudy_trend;
vdudy_trend2(vdudy_sig==0)=NaN;

subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),vdudy_trend2,100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-1e-8 1e-8]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('trend in v du/dy (m s^{-1}y^{-1})')
axis([140 280 1 size(depth_iso,3)])


delta_iso_vdudy = squeeze(nanmean(nanmean(vdudy(:,I_lat,:,end-29:end),4),2))- squeeze(nanmean(nanmean(vdudy(:,I_lat,:,1:30),4),2));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_vdudy',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta v du/dy (m s^{-2})')
axis([140 280 1 size(depth_iso,3)])


%%
close all

wdudz_trend2 = wdudz_trend;
wdudz_trend2(wdudz_sig==0)=NaN;

figure(4)
subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),wdudz_trend2,100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-1e-8 1e-8]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('trend in w du/dz (m s^{-1}y^{-1})')
axis([140 280 1 size(depth_iso,3)])

delta_iso_wdudz = squeeze(nanmean(nanmean(wdudz(:,I_lat,:,end-29:end),4),2))- squeeze(nanmean(nanmean(wdudz(:,I_lat,:,1:30),4),2));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_wdudz',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta w du/dz (m s^{-2})')
axis([140 280 1 size(depth_iso,3)])

%%
close all


corf_trend2 = corf_trend;
corf_trend2(corf_sig==0)=NaN;

figure(5)
subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),corf_trend2,100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-1e-8 1e-8]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('trend in 2 \Omega sin(\phi) (m s^{-1}y^{-1}))')
axis([140 280 1 size(depth_iso,3)])


delta_iso_corf = squeeze(nanmean(nanmean(corf(:,I_lat,:,end-29:end),4),2))- squeeze(nanmean(nanmean(corf(:,I_lat,:,1:30),4),2));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_corf',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta 2 \Omega sin(\phi) (m s^{-2})')
axis([140 280 1 size(depth_iso,3)])

%%

close all

frich_trend2 = frich_trend;
frich_trend2(frich_sig==0)=NaN;

figure(6)
subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),frich_trend2,100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-1e-8 1e-8]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('trend in A_h (d^2u/dx^2 + d^2u/dy^2) (m s^{-1}y^{-1}))')
axis([140 280 1 size(depth_iso,3)])

delta_iso_frich = squeeze(nanmean(nanmean(frich(:,I_lat,:,end-29:end),4),2))- squeeze(nanmean(nanmean(frich(:,I_lat,:,1:30),4),2));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_frich',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta A_h (d^2u/dx^2 + d^2u/dy^2) (m s^{-2})')
axis([140 280 1 size(depth_iso,3)])


%%

close all

fricv_trend2 = fricv_trend;
fricv_trend2(fricv_sig==0)=NaN;

figure(7)
subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),fricv_trend2,100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-1e-8 1e-8]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('trend in d/dz(A_v*du/dz) (m s^{-1}y^{-1}))')
axis([140 280 1 size(depth_iso,3)])

delta_iso_fricv = squeeze(nanmean(nanmean(fricv(:,I_lat,:,end-29:end),4),2))- squeeze(nanmean(nanmean(fricv(:,I_lat,:,1:30),4),2));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_fricv',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta d/dz(A_v*du/dz) (m s^{-2})')
axis([140 280 1 size(depth_iso,3)])

%%
close all

zpgf_trend2 = zpgf_trend;
zpgf_trend2(zpgf_sig==0)=NaN;

figure(7)
subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),zpgf_trend2,100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-1e-8 1e-8]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('trend in -\rho^{-1} dp/dx (m s^{-2}) (m s^{-1}y^{-1}))')
axis([140 280 1 size(depth_iso,3)])

delta_iso_zpgf = squeeze(nanmean(nanmean(zpgf(:,I_lat,:,end-29:end),4),2))- squeeze(nanmean(nanmean(zpgf(:,I_lat,:,1:30),4),2));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_zpgf',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth_iso,3),squeeze(nanmean(nanmean(u_iso(:,I_lat,:,:),4),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta -\rho^{-1} dp/dx (m s^{-2}) (m s^{-2})')
axis([140 280 1 size(depth_iso,3)])