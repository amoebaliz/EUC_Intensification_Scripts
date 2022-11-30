close all
clear all
clc

load('ashort.mat')

% RED/BLUE COLORMAP DEFINITION
colorStart = [0 0 1];           % Blue (negative end of spectrum)
colorCenter = [1 1 1];    % Grey (to be centered at 'caxis' 0)
colorEnd = [1 0 0];             % Red  (positive end of spectrum)
num = 30;

cmap1 = zeros(num,3);
cmap2 = cmap1;

for j = 1:3
    
    cmap1(1:num,j)= linspace(colorStart(j), colorCenter(j),num);
    cmap2(1:num,j) = linspace(colorCenter(j), colorEnd(j),num);
end
cmap = [cmap1(1:end-1,:);cmap2(:,:)];

depth_short= depth(1:size(u_short,3));

figure(1)
colormap(cmap)
[~,h]=contourf(slon,depth(1:14),squeeze(mean(mean(u_short(:,2:3,:,1),2),4))',500);
set(h,'EdgeColor','none')
set(gca,'YDir','reverse')

nlon = size(u_short,1);
nlat = size(u_short,2);
ndep = size(u_short,3);


u_trend=NaN(nlon,nlat,ndep); 
u_sig = NaN(nlon,nlat,ndep); 

for j = 1:nlon
    for k = 1:nlat
        for d = 1:ndep
            [u_trend(j,k,d),~,u_sig(j,k,d)] = trend_stat(squeeze(u_short(j,k,d,:)),95);
        end
    end
end

u_trend(u_sig==0)=NaN;

figure(2)
colormap(cmap)
[~,h]=contourf(slon,depth(1:14),squeeze(mean(u_trend(:,2:3,:),2))',500);
set(h,'EdgeColor','none')
set(gca,'YDir','reverse')

%%

% define constants

dx=111320/2;
dy=110574/2;
dz=2;
rho0=1028;
g=-9.81;
omega=7.292*10^(-5);


% compute density
rho=NaN(size(temp_short));
for z=1:ndep
    rho(:,:,z)=sw_dens(squeeze(salt_short(:,:,z)),squeeze(temp_short(:,:,z)),depth_short(z));
end; clear z

% compute pressure
p=NaN(size(temp_short));
for x=1:nlon
    for y=1:nlat
        for z=1:ndep
            p(x,y,z)=-1*g*(sum(rho(x,y,1:z),3)*dz+ssh_short(x,y)*rho(x,y,1));
        end
    end
end; clear x y z rho

% compute zonal pressure gradient force

zpgf=NaN(size(u_short));
dpdx=NaN(size(u_short));
for x=2:nlon-1
    for y=1:nlat
        for z=1:ndep
            dpdx(x,y,z)=(p(x+1,y,z)-p(x-1,y,z))/(2*dx);
            zpgf(x,y,z)=(-1/rho0)*dpdx(x,y,z);
        end
    end
end; clear x y z dpdx


% compute nonlinear advective terms

ududx=NaN(size(u_short));
vdudy=NaN(size(u_short));
wdudz=NaN(size(u_short));

for x=2:nlon-1
    for y=1:nlat
        for z=1:ndep
            ududx(x,y,z)=u_short(x,y,z)*(u_short(x+1,y,z)-u_short(x-1,y,z))/(2*dx);
        end
    end
end; clear x y z

for x=1:nlon
    for y=2:nlat-1
        for z=1:ndep
            vdudy(x,y,z)=v_short(x,y,z)*(u_short(x,y+1,z)-u_short(x,y-1,z))/(2*dy);
        end
    end
end; clear x y z

for x=1:nlon
    for y=1:nlat
        for z=2:ndep-1
            wdudz(x,y,z)=w_short(x,y,z)*(u_short(x,y,z-1)-u_short(x,y,z+1))/(2*dz);
        end
    end
end; clear x y z

figure(3)
colormap(cmap)
[~,h]=contourf(slon(2:end-1),depth(2:13),squeeze(nanmean(nanmean(ududx(2:end-1,2:3,2:13,:),2),4))',500);
set(h,'EdgeColor','none')
set(gca,'YDir','reverse')
colorbar

% center diff u

center_diff=(u_short(:,:,:,3:end)-u_short(:,:,:,1:end-2))/2;

% convert to m/s second-1
accel = center_diff/(30*24*60*60);

figure(4)
colormap(cmap)
[~,h]=contourf(slon(2:end-1),depth(2:13),squeeze(nanmean(nanmean(accel(2:end-1,2:3,2:13,:),2),4))',500);
set(h,'EdgeColor','none')
set(gca,'YDir','reverse')
colorbar



% compute coriolis force

corf=NaN(size(u_short));

for x=1:nlon
    for y=1:nlat
        for z=1:ndep
            corf(x,y,z)=2*omega*sind(slat(y))*v_short(x,y,z);
        end
    end
end; clear x y z


% compute potential density

rhoth=NaN(size(temp_short));
for z=1:ndep
        rhoth(:,:,z)=sw_pden(squeeze(salt_short(:,:,z)),squeeze(temp_short(:,:,z)),depth(z),0);
end; clear z salt


% VERTICALLY VARYING COEFFICIENT OF VERTICAL EDDY VISCOSITY (Av) SEE QIAO
% AND WEISBERG (1997, JPO, FIG. 14) THERE ARE TWO OPTIONS: TIE INFLECTION
% POINT TO THE THERMOCLINE (MAX DT/DZ) OR TO THE EUC (MAX U)

Av=NaN(size(u));
Av_sfc=45*10^(-4); Av_tcline=3*10^(-4); Av_deep=15*10^(-4);





dudt = zpgf+ududx+vdudy+wdudz;

disp('max')

max(max(max(max(accel))))
max(max(max(max(dudt))))

max(max(max(max(zpgf))))
max(max(max(max(ududx))))
max(max(max(max(vdudy))))
max(max(max(max(wdudz))))

disp('min')
min(min(min(min(accel))))
min(min(min(min(dudt))))

min(min(min(min(zpgf))))
min(min(min(min(ududx))))
min(min(min(min(vdudy))))
min(min(min(min(wdudz))))
