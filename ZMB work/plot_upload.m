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

save june_isopyc ...
    ...
    lon_iso lat_iso depth_iso u_iso dudt ududx vdudy wdudz zpgf corf


% plotting the trends and the average u_iso for each of the saved fields

I_lat = ceil(length(lat_iso)/2);

% allocate space in advance
u_iso_trend = NaN(length([5:dz:400]'),length(lon_iso));

u_iso_plusminus = u_iso_trend; u_iso_sig = u_iso_trend; 
ududx_trend = u_iso_trend; ududx_plusminus = u_iso_trend; ududx_sig = u_iso_trend;
vdudy_trend = u_iso_trend; vdudy_plusminus = u_iso_trend; vdudy_sig = u_iso_trend;
wdudz_trend = u_iso_trend; wdudz_plusminus = u_iso_trend; wdudz_sig = u_iso_trend;

zpgf_trend = u_iso_trend; zpgf_plusminus = u_iso_trend; zpgf_sig = u_iso_trend;
corf_trend = u_iso_trend; corf_plusminus = u_iso_trend; corf_sig = u_iso_trend;

fricv_trend = u_iso_trend; fricv_plusminus = u_iso_trend; fricv_sig = u_iso_trend;
frich_trend = u_iso_trend; frich_plusminus = u_iso_trend; frich_sig = u_iso_trend;

for j = 1:length(lon_iso)
    for k = 1:length([5:dz:400]')
        
       [u_iso_trend(k,j),u_iso_plusminus(k,j),u_iso_sig(k,j)] = trend(squeeze(u_iso(j,I_lat,k,:)),99);
       
       [ududx_trend(k,j),ududx_plusminus(k,j),ududx_sig(k,j)] = trend(squeeze(ududx(j,I_lat,k,:)),99);
       [vdudy_trend(k,j),vdudy_plusminus(k,j),vdudy_sig(k,j)] = trend(squeeze(vdudy(j,I_lat,k,:)),99);
       [wdudz_trend(k,j),wdudz_plusminus(k,j),wdudz_sig(k,j)] = trend(squeeze(wdudz(j,I_lat,k,:)),99);
       
       [zpgf_trend(k,j),zpgf_plusminus(k,j),zpgf_sig(k,j)] = trend(squeeze(zpgf(j,I_lat,k,:)),99);
       [corf_trend(k,j),corf_plusminus(k,j),corf_sig(k,j)] = trend(squeeze(corf(j,I_lat,k,:)),99);
       
       [fricv_trend(k,j),fricv_plusminus(k,j),fricv_sig(k,j)] = trend(squeeze(fricv(j,I_lat,k,:)),99);
       [frich_trend(k,j),frich_plusminus(k,j),frich_sig(k,j)] = trend(squeeze(frich(j,I_lat,k,:)),99);
        
    end
end

figure(1)

subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),u_iso_trend',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('u du/dx (m s^{-2})')
axis([140 280 1 size(depth,3)])

delta_iso_u = squeeze(mean(iso_u(:,I_lat,:,end-29:end)))- squeeze(mean(iso_u(:,I_lat,:,1:30)));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_u',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta u du/dx (m s^{-2})')
axis([140 280 1 size(depth,3)])



figure(2)
subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),ududx_trend',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('v du/dy (m s^{-2})')
axis([140 280 1 size(depth,3)])


delta_iso_ududx = squeeze(mean(ududx(:,I_lat,:,end-29:end)))- squeeze(mean(ududx(:,I_lat,:,1:30)));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_ududx',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta v du/dy (m s^{-2})')
axis([140 280 1 size(depth,3)])







figure(3)

subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),vdudy_trend',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('w du/dz (m s^{-2})')
axis([140 280 1 size(depth,3)])


delta_iso_vdudy = squeeze(mean(vdudy(:,I_lat,:,end-29:end)))- squeeze(mean(vdudy(:,I_lat,:,1:30)));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_vdudy',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta w du/dz (m s^{-2})')
axis([140 280 1 size(depth,3)])







figure(4)
subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),wdudz_trend',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('-\rho^{-1} dp/dx (m s^{-2})')
axis([140 280 1 size(depth,3)])

delta_iso_wdudz = squeeze(mean(wdudz(:,I_lat,:,end-29:end)))- squeeze(mean(wdudz(:,I_lat,:,1:30)));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_wdudz',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta -\rho^{-1} dp/dx (m s^{-2})')
axis([140 280 1 size(depth,3)])




figure(5)
subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),corf_trend',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('2 \Omega sin(\phi) (m s^{-2})')
axis([140 280 1 size(depth,3)])


delta_iso_corf = squeeze(mean(corf(:,I_lat,:,end-29:end)))- squeeze(mean(corf(:,I_lat,:,1:30)));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_corf',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta 2 \Omega sin(\phi) (m s^{-2})')
axis([140 280 1 size(depth,3)])


figure(6)
subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),frich_trend',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('A_h (d^2u/dx^2 + d^2u/dy^2) (m s^{-2})')
axis([140 280 1 size(depth,3)])

delta_iso_frich = squeeze(mean(frich(:,I_lat,:,end-29:end)))- squeeze(mean(frich(:,I_lat,:,1:30)));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_frich',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta A_h (d^2u/dx^2 + d^2u/dy^2) (m s^{-2})')
axis([140 280 1 size(depth,3)])




figure(7)
subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),fricv_trend',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('d/dz(A_v*du/dz) (m s^{-2})')
axis([140 280 1 size(depth,3)])

delta_iso_fricv = squeeze(mean(fricv(:,I_lat,:,end-29:end)))- squeeze(mean(fricv(:,I_lat,:,1:30)));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_fricv',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta d/dz(A_v*du/dz) (m s^{-2})')
axis([140 280 1 size(depth,3)])


figure(8)
subplot(2,1,2)
[~,h]=contourf(lon,1:size(depth_iso,3),squeeze(mean(dudt(:,lat>=-0.5&lat<=0.5,:),2))',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('du/dt (sum) (m s^{-2})')
axis([140 280 1 size(depth,3)])


delta_iso_dudt = squeeze(mean(dudt(:,I_lat,:,end-29:end)))- squeeze(mean(dudt(:,I_lat,:,1:30)));

subplot(2,1,1)
[~,h]=contourf(lon,1:size(depth_iso,3),delta_iso_fricv',100);
set(h,'EdgeColor','none')
colorbar; colormap(redblue); caxis([-5e-7 5e-7]);
set(gca,'YDir','reverse')
hold on; contour(lon,1:size(depth,3),squeeze(mean(u_iso(:,lat>=-0.5&lat<=0.5,:),2))',[0.5 0.5],'k','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Longitude'); ylabel('Layer'); title('\Delta du/dt (sum) (m s^{-2})')
axis([140 280 1 size(depth,3)])
