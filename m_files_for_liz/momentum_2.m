close all
clear all
%% Load relevant stuff
load momentum.mat
load filled_contours.mat
%% initialize terms

[a b c] = size(avg_u);

udu_term = zeros(a-2,b-2,c-2);
vdu_term = udu_term;
wdu_term = udu_term;
dudz = zeros([1,c-2]);

I_lat =  find(lat<=0.5 & lat>=-0.5);

%% calculate all terms 

% for h = 1:12
   for j = 2:(length(lon)-1)
       for k = 2:(length(lat)-1)
            for d = 2:(length(depth)-1)
           
%% Calculate udu/dx
        u_diff_x = avg_u(j+1,k,d)-avg_u(j-1,k,d);                   % calculat delta u
        udu_term(j-1,k-1,d-1) = avg_u(j,k,d)*(u_diff_x/delta_Lon);  % u*delt_u/delta_x

%% Calculate vdu/dy
        u_diff_y = avg_u(j,k+1,d)-avg_u(j,k-1,d);
        vdu_term(j-1,k-1,d-1) = avg_v(j,k,d)*(u_diff_y/delta_Lat);
%% Calculate wdu/dz
        u_diff_z = avg_u(j,k,d-1)-avg_u(j,k,d+1);
        wdu_term(j-1,k-1,d-1) = avg_w(j,k,d)*(u_diff_z/delta_z_2(d-1));

%% Calculate Fx (windstress and propagated to depth)
        % get du/dz for all depths before ending depth for loop
        dudz(d-1) = u_diff_z/delta_z_2(d-1);
          
           end
       for dd = 2:17
            Fx(j-1,k-1,dd-1)=(1/rho)*mu*((dudz(dd-1)-dudz(dd+1))/delta_z_3(dd-1));
       end
%% Calculate dP/dx
        
% the horizontal pressure gradient is identical everywhere within the
% fluid; so it doesn't matter how deep you're looking

% this simplifies to: g*deltaZ/deltaX
% g = gravitational constant (9.81 m s^-2)
% deltaX (distance in longitude is 1 degree and assumed constant: 111.320 km)
% deltaZ (difference in sea surface height) requires SSH at lon(i-1), (i+1)    
        delta_SSH = avg_ssh(j+1,k)-avg_ssh(j-1,k); 
        pressure_term(j-1,k-1) = -1*grav*delta_SSH/delta_Lon;
       end
    end
%%
% close all
 cmin = -.0000003;
 cmax =  .0000003;
% 
 fig1 = figure;
 set(fig1,'Position',[100 100 1200 800]); 

% for k = 1:12
%     
%     
 subplot(5,1,1)
 m = pcolor(repmat(lon(2:319)',[18,1]),repmat(depth(2:19),[1,318]),double(repmat((mean(pressure_term(:,I_lat),2))',[18,1])));
 set(gca,'YDir','reverse')
 set(m, 'EdgeColor', 'none')
 caxis([cmin cmax])
 title('Pressure Gradient Force')    
%     
 subplot(5,1,2)
 m = pcolor(repmat(lon(2:319)',[16,1]),repmat(depth(3:18),[1,318]),double(squeeze(nanmean(Fx(:,I_lat,:),2)))');
 set(gca,'YDir','reverse')
 set(m, 'EdgeColor', 'none')
 caxis([cmin cmax])
 title('Fx')    
% 
 subplot(5,1,3)
 m = pcolor(repmat(lon(2:319)',[18,1]),repmat(depth(2:19),[1,318]),double(squeeze(nanmean(udu_term(:,I_lat,:),2)))');
 set(gca,'YDir','reverse')
 set(m, 'EdgeColor', 'none')
 caxis([cmin cmax])
title('ududx')

 subplot(5,1,4)
 p = pcolor(repmat(lon(2:319)',[18,1]),repmat(depth(2:19),[1,318]),double(squeeze(nanmean(vdu_term(:,I_lat,:),2)))');
% %contourf(lon(21:319),depth(2:19),wdu_term(20:end,:,k)')
 set(gca,'YDir','reverse')
% %cBar = colorbar;
 set(p, 'EdgeColor', 'none')
 caxis([cmin cmax])
title('vdudy')
% 
 subplot(5,1,5)
 t = pcolor(repmat(lon(2:319)',[18,1]),repmat(depth(2:19),[1,318]),double(squeeze(nanmean(wdu_term(:,I_lat,:),2)))');
% %contourf(lon(21:319),depth(2:19),wdu_term(20:end,:,k)')
 set(gca,'YDir','reverse')
% %cBar = colorbar;
caxis([cmin cmax])
 set(t, 'EdgeColor', 'none')
 title('wdudz')

% 
%%
P_term = double(repmat((mean(pressure_term(:,I_lat),2))',[18,1]));
Fx_term = double(squeeze(nanmean(Fx(:,I_lat,:),2)))';

u_term = double(squeeze(nanmean(udu_term(:,I_lat,:),2)))';
v_term = double(squeeze(nanmean(vdu_term(:,I_lat,:),2)))';
w_term = double(squeeze(nanmean(wdu_term(:,I_lat,:),2)))';

sum_term = (P_term(2:17,:)+ Fx_term) - (u_term(2:17,:)+v_term(2:17,:)+w_term(2:17,:));

T = double(squeeze(nanmean(avg_u(:,I_lat,:),2)))';
T = T*100;                                            % centimeters/second

sum_term(isnan(euc_trend)) = NaN;

%COLORMAP DEFINITION
colorStart = [0 .1 .9];           % Blue (negative end of spectrum)
colorCenter = [1 1 1];    % Grey (to be centered at 'caxis' 0)
colorEnd = [.9 0 0];             % Red  (positive end of spectrum)
num = 30;

cmap1 = zeros(num,3);
cmap2 = cmap1;

for j = 1:3
    
    cmap1(1:num,j)= linspace(colorStart(j), colorCenter(j),num);
    cmap2(1:num,j) = linspace(colorCenter(j), colorEnd(j),num);
end

fig1 = figure('Color',[1 1 1]);

set(fig1,'Position',[400 600 1100 600]) 
cmap = [cmap1(1:end-1,:);cmap2(:,:)];
colormap(cmap)

t = pcolor(repmat(lon(2:319)',[16,1]),repmat(depth(3:18),[1,318]),double(sum_term));
%t = contourf(lon(2:319),depth(3:18),double(sum_term));

hold on
qq = contour(lon(2:319),depth(3:18),T(3:18,2:319),[-30,-15],'-.b','LineWidth',2); 
q = contour(lon(2:319),depth(3:18),T(3:18,2:319),[0, 0, 0],':k', 'LineWidth',1);
qqq = contour(lon(2:319),depth(3:18),T(3:18,2:319),[15,30,45,60],'k', 'LineWidth',2);
set(gca,'YDir','reverse')
set(t, 'EdgeColor', 'none')
caxis([-.0000008 .0000008])
xlim([150 270])
ylim([depth(3) depth(18)])
%%
figure(3)
t = pcolor(repmat(lon(2:319)',[16,1]),repmat(depth(3:18),[1,318]),double(squeeze(nanmean(Fx(:,9:10,:),2)))');
 caxis([-.0000000001 .0000000001])
set(gca,'YDir','reverse')
set(t, 'EdgeColor', 'none')
xlim([150 270])

%%
% 
% pause(2)
% end
% %%
% 
% close all
% fig2 = figure;
% L = pcolor(repmat(lon(2:319)',[12,1]),repmat((1:12)',[1,318]),pressure_term);
% set(L, 'EdgeColor', 'none')
% cmin = -0.0000006;
% cmax =  0.0000006;
% caxis([cmin cmax])