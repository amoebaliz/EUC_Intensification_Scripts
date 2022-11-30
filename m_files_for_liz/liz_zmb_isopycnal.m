function [lon,lat,depth,u,dudt,ududx,vdudy,wdudz,zpgf,corf,frich,fricv] = liz_zmb_isopycnal(lon,lat,depth,u,v,w,temp,salt,ssh)

% Isopycnal zonal momentum budget with spatially varying vertical eddy
% viscosity for the equatorial Pacific Ocean in extended SODA

%% define constants

dx=111320/2;
dy=110574/2;
dz=2;
rho0=1028;
g=-9.81;
omega=7.292*10^(-5);
Ah=1.5*10^3; % Wallcroft et al. (2005, GRL)
% For Av, see constant vs. varying options below

%% interpolate to uniform depth grid

depth(1)=5.0;
depth2=(5:dz:400)';

u2 = permute(interp1(depth,permute(u,[3 1 2 4]),depth2),[2 3 1 4]);
v2 = permute(interp1(depth,permute(v,[3 1 2 4]),depth2),[2 3 1 4]);
w2 = permute(interp1(depth,permute(w,[3 1 2 4]),depth2),[2 3 1 4]);
temp2 = permute(interp1(depth,permute(temp,[3 1 2 4]),depth2),[2 3 1 4]);
salt2 = permute(interp1(depth,permute(salt,[3 1 2 4]),depth2),[2 3 1 4]);

warning('on','all');

u=u2; clear u2
v=v2; clear v2
w=w2; clear w2
temp=temp2; clear temp2
salt=salt2; clear salt2
depth=depth2; clear depth2

%% interpolate to lat grid with an equatorial gridpoint

lat2=(-4.5:0.5:4.5)';

warning('off','all');

u2 = permute(interp1(lat,permute(u,[2 1 3 4]),lat2),[2 1 3 4]);
v2 = permute(interp1(lat,permute(v,[2 1 3 4]),lat2),[2 1 3 4]);
w2 = permute(interp1(lat,permute(w,[2 1 3 4]),lat2),[2 1 3 4]);
temp2 = permute(interp1(lat,permute(temp,[2 1 3 4]),lat2),[2 1 3 4]);
salt2 = permute(interp1(lat,permute(salt,[2 1 3 4]),lat2),[2 1 3 4]);
ssh2 = permute(interp1(lat,permute(ssh,[2 1 3]),lat2),[2 1 3]);

warning('on','all');

u=u2; clear u2
v=v2; clear v2
w=w2; clear w2
temp=temp2; clear temp2
salt=salt2; clear salt2
ssh=ssh2; clear ssh2
lat=lat2; clear lat2

%% compute nonlinear advective terms

ududx = u(2:end-1,:,:,:).*(u(3:end,:,:,:)-u(1:end-2,:,:,:))/(2*dx);
vdudy = v(:,2:end-1,:,:).*(u(:,3:end,:,:)-u(:,1:end-2,:,:))/(2*dy);
wdudz = w(:,:,2:end-1,:).*(u(:,:,3:end,:)-u(:,:,1:end-2,:))/(2*dz);

%% compute coriolis force

lat_inf = repmat(sind(lat)',[320 1 20 1716]);
corf = 2*omega*lat_inf.*v;

%% time itterated calculations
zpgf=NaN(size(u));

for t=1:length(time)

% compute density

rho=NaN(size(temp));
for z=1:length(depth)
    rho(:,:,z)=sw_dens(squeeze(salt(:,:,z,t)),squeeze(temp(:,:,z,t)),depth(z));
end; clear z

% compute pressure

p=NaN(size(temp));
for x=1:length(lon)
    for y=1:length(lat)
        for z=1:length(depth)
            p(x,y,z)=-1*g*(sum(rho(x,y,1:z),3)*dz+ssh(x,y,t)*rho(x,y,1));
        end
    end
end; clear x y z rho

% compute zonal pressure gradient force

dpdx=NaN(size(u,1:3));
for x=2:length(lon)-1
    for y=1:length(lat)
        for z=1:length(depth)
            dpdx(x,y,z)=(p(x+1,y,z)-p(x-1,y,z))/(2*dx);
            zpgf(x,y,z,t)=(-1/rho0)*dpdx(x,y,z);
        end
    end
end; clear x y z dpdx

% compute potential density

rhoth=NaN(size(temp));
for z=1:length(depth)
        rhoth(:,:,z)=sw_pden(squeeze(salt(:,:,z)),squeeze(temp(:,:,z)),depth(z),0);
end; clear z salt


%%


% Constant Av options

% Av=ones(size(u))*1.66*10^(-3); % Bryden & Brady (1985)
% Av=ones(size(u))*1*10^(-4); % MITgcm manual

% VERTICALLY VARYING COEFFICIENT OF VERTICAL EDDY VISCOSITY (Av) SEE QIAO
% AND WEISBERG (1997, JPO, FIG. 14) THERE ARE TWO OPTIONS: TIE INFLECTION
% POINT TO THE THERMOCLINE (MAX DT/DZ) OR TO THE EUC (MAX U)

Av=NaN(size(u));
Av_sfc=45*10^(-4); Av_tcline=3*10^(-4); Av_deep=15*10^(-4);

% OPTION 1: TIE INFLECTION POINT TO THE THERMOCLINE

% deep is depth of thermocline + 100 m. note that QW97 says the values
% below the EUC might be at least an order of magnitude smaller than this
% based on measurements.

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
end; clear x y z_tcline z_deep a b c d Av_sfc Av_tcline Av_deep temp


% convert to isopycnal layers

isopyc_layer_bnds=1020:0.2:1028;

% non-uniform layers
% nlayers=40;
% isopyc_layer_bnds=zeros(nlayers+1,1);
% isopyc_layer_bnds(1)=1020;
% for layer=2:nlayers+1
%     isopyc_layer_bnds(layer)=isopyc_layer_bnds(layer-1)+layer/107;
% end; clear nlayers layer
% isopyc_layer_bnds

depth_iso=NaN(length(lon),length(lat),length(isopyc_layer_bnds)-1);
Av_iso=NaN(length(lon),length(lat),length(isopyc_layer_bnds)-1);
p_iso=NaN(length(lon),length(lat),length(isopyc_layer_bnds)-1);
u_iso=NaN(length(lon),length(lat),length(isopyc_layer_bnds)-1);
v_iso=NaN(length(lon),length(lat),length(isopyc_layer_bnds)-1);
w_iso=NaN(length(lon),length(lat),length(isopyc_layer_bnds)-1);
rhoth_iso=NaN(length(lon),length(lat),length(isopyc_layer_bnds)-1);
zpgf_iso=NaN(length(lon),length(lat),length(isopyc_layer_bnds)-1);
ududx_iso=NaN(length(lon),length(lat),length(isopyc_layer_bnds)-1);
vdudy_iso=NaN(length(lon),length(lat),length(isopyc_layer_bnds)-1);
wdudz_iso=NaN(length(lon),length(lat),length(isopyc_layer_bnds)-1);
corf_iso=NaN(length(lon),length(lat),length(isopyc_layer_bnds)-1);

for x=1:length(lon)
    for y=1:length(lat)
        for layer=1:length(isopyc_layer_bnds)-1
            zrange=depth(squeeze(rhoth(x,y,:))>=isopyc_layer_bnds(layer)&squeeze(rhoth(x,y,:))<=isopyc_layer_bnds(layer+1));
            if isfinite(zrange)
                depth_iso(x,y,layer)=nanmean(depth(depth>=zrange(1)&depth<=zrange(length(zrange))));
                Av_iso(x,y,layer)=nanmean(Av(x,y,depth>=zrange(1)&depth<=zrange(length(zrange))),3);
                p_iso(x,y,layer)=nanmean(p(x,y,depth>=zrange(1)&depth<=zrange(length(zrange))),3);
                u_iso(x,y,layer)=nanmean(u(x,y,depth>=zrange(1)&depth<=zrange(length(zrange))),3);
                v_iso(x,y,layer)=nanmean(v(x,y,depth>=zrange(1)&depth<=zrange(length(zrange))),3);
                w_iso(x,y,layer)=nanmean(w(x,y,depth>=zrange(1)&depth<=zrange(length(zrange))),3);
                rhoth_iso(x,y,layer)=nanmean(rhoth(x,y,depth>=zrange(1)&depth<=zrange(length(zrange))),3);
                zpgf_iso(x,y,layer)=nanmean(zpgf(x,y,depth>=zrange(1)&depth<=zrange(length(zrange))),3);
                ududx_iso(x,y,layer)=nanmean(ududx(x,y,depth>=zrange(1)&depth<=zrange(length(zrange))),3);
                vdudy_iso(x,y,layer)=nanmean(vdudy(x,y,depth>=zrange(1)&depth<=zrange(length(zrange))),3);
                wdudz_iso(x,y,layer)=nanmean(wdudz(x,y,depth>=zrange(1)&depth<=zrange(length(zrange))),3);
                corf_iso(x,y,layer)=nanmean(corf(x,y,depth>=zrange(1)&depth<=zrange(length(zrange))),3);
            end
        end
    end
end; clear x y layer zrange

% % check conversion
% 
% figure(100)
% subplot(2,2,1)
% [~,h]=contourf(lon,depth,squeeze(rhoth(:,lat==0,:))',50);
% set(h,'EdgeColor','none')
% set(gca,'YDir','reverse')
% colorbar
% hold on
% contour(lon,depth,squeeze(rhoth(:,lat==0,:))',isopyc_layer_bnds,'k','LineWidth',1)
% subplot(2,2,2)
% contourf(lon,depth,squeeze(u(:,lat==0,:))',[-1:0.2:1])
% set(gca,'YDir','reverse')
% colorbar
% subplot(2,2,3)
% [~,h]=contourf(lon,1:length(isopyc_layer_bnds)-1,squeeze(rhoth_iso(:,lat==0,:))',50);
% set(h,'EdgeColor','none')
% set(gca,'YDir','reverse')
% colorbar
% hold on
% contour(lon,1:length(isopyc_layer_bnds)-1,squeeze(rhoth_iso(:,lat==0,:))',isopyc_layer_bnds,'k','LineWidth',1)
% subplot(2,2,4)
% contourf(lon,1:length(isopyc_layer_bnds)-1,squeeze(u_iso(:,lat==0,:))',[-1:0.2:1])
% set(gca,'YDir','reverse')
% colorbar

clear rhoth rhoth_iso

depth=depth_iso; clear depth_iso
Av=Av_iso; clear Av_iso
p=p_iso; clear p_iso
u=u_iso; clear u_iso
v=v_iso; clear v_iso
w=w_iso; clear w_iso
zpgf=zpgf_iso; clear zpgf_iso
ududx=ududx_iso; clear ududx_iso
vdudy=vdudy_iso; clear vdudy_iso
wdudz=wdudz_iso; clear wdudz_iso
corf=corf_iso; clear corf_iso





% compute horizontal friction terms

fricx=NaN(size(u));
fricy=NaN(size(u));

dudx=NaN(size(u));
for x=2:length(lon)-1
    for y=1:length(lat)
        for z=1:size(depth,3)
            dudx(x,y,z)=(u(x+1,y,z)-u(x-1,y,z))/(2*dx);
        end
    end
end; clear x y z

d2udx2=NaN(size(u));
for x=3:length(lon)-2
    for y=1:length(lat)
        for z=1:size(depth,3)
            d2udx2(x,y,z)=(dudx(x+1,y,z)-dudx(x-1,y,z))/(2*dx);
            fricx(x,y,z)=Ah*d2udx2(x,y,z);
        end
    end
end; clear x y z

dudy=NaN(size(u));
for x=1:length(lon)
    for y=2:length(lat)-1
        for z=1:size(depth,3)
            dudy(x,y,z)=(u(x,y+1,z)-u(x,y-1,z))/(2*dy);
        end
    end
end; clear x y z

d2udy2=NaN(size(u));
for x=1:length(lon)
    for y=3:length(lat)-2
        for z=1:size(depth,3)
            d2udy2(x,y,z)=(dudy(x,y+1,z)-dudy(x,y-1,z))/(2*dy);
            fricy(x,y,z)=Ah*d2udy2(x,y,z);
        end
    end
end; clear x y z

% sum fricx and fricy to get frich

frich=fricx+fricy;





% compute vertical friction term

fricv=NaN(size(u));

dudz=NaN(size(u));
Av_dudz=NaN(size(u));
for x=1:length(lon)
    for y=1:length(lat)
        for z=2:size(depth,3)-1
            dudz(x,y,z)=(u(x,y,z-1)-u(x,y,z+1))/( depth(x,y,z+1)-depth(x,y,z-1) );
            Av_dudz(x,y,z)=Av(x,y,z)*dudz(x,y,z);
        end
    end
end; clear x y z

for x=1:length(lon)
    for y=1:length(lat)
        for z=3:size(depth,3)-2
            fricv(x,y,z)=(Av_dudz(x,y,z-1)-Av_dudz(x,y,z+1))/( depth(x,y,z+1)-depth(x,y,z-1) );
        end
    end
end; clear x y z


end


%% compute the sum, i.e., dudt

dudt = -1*ududx + -1*vdudy + -1*wdudz + zpgf + corf + frich + fricv;

