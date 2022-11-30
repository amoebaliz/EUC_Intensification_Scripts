close all
clear all
clc

lon = ncread('SODA_2.2.6_Trop_Pac_WE.cdf','LON241_580');
lat = ncread('SODA_2.2.6_Trop_Pac_WE.cdf','LAT142_161');
time = ncread('SODA_2.2.6_Trop_Pac_u.cdf','TIME');
depth = ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');

we = ncread('SODA_2.2.6_Trop_Pac_WE.cdf','WE');
taux = ncread('SODA_2.2.6_Trop_Pac_taux.cdf','TAUX_ENS_MN');
u = ncread('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN');

% removing missing values -1e+34
we(we<-10000) = NaN;
taux(taux<-10000) = NaN;
u(u<-10000) = NaN;

% Indexing spatial bounds
 I_lat = find(lat<=0.5 & lat>=-0.5);     
 I_lon = find(lon>=149 &lon<=271);
 I_time = find(time>1 & time<time(end));
 
 slat = lat(I_lat);
 slon = lon(I_lon);
 stime = time(I_time);
 
 nt = length(time);
 nlon = length(slon);
 nlat= length(slat);
 ndep = length(depth);
 
 %%
 
 we_avg = squeeze(mean(mean(we(I_lon, I_lat,I_time),3),2));
 we_time_mean = squeeze(mean(mean(we(I_lon, I_lat,I_time),2),1));
 taux_avg = squeeze(mean(mean(taux(I_lon, I_lat,I_time),3),2));
 u_avg = squeeze(mean(mean(u(I_lon, I_lat,:,I_time),4),2));
 
 %%
 subplot(2,1,1)
 plot(slon, we_avg)
 
 subplot(2,1,2)
 plot(stime, we_time_mean)
 
 % sum(isnan(we_time_mean))
 
 %%

 we_trend = NaN(nlon,1); we_plusminus = we_trend; we_sig = we_trend;
 taux_trend = NaN(nlon,1); taux_plusminus = taux_trend; taux_sig = taux_trend;
 u_trend = NaN(nlon, ndep); u_sig = u_trend;
 
 for j = 1:nlon
    
     [we_trend(j), we_plusminus(j), we_sig(j)] =trend_stat(squeeze(mean(we(I_lon(j),I_lat,I_time),2)),95);
     [taux_trend(j), taux_plusminus(j), taux_sig(j)] =trend_stat(squeeze(mean(taux(I_lon(j),I_lat,I_time),2)),95);
    
     for k = 1:ndep
         
         [u_trend(j,k), ~ , u_sig(j,k)] = trend_stat(squeeze(mean(u(I_lon(j),I_lat,k,I_time),2)),95);
         
     end
 
 end
 
 % we_trend(we_sig==0)=NaN;
 %%
 close all
 
 figure('Position', [300 10 800 1000])
set(gcf, 'Renderer', 'Painters')

load redblue
colormap(redblue)
 
hsub = subplot(3,1,1);
set(hsub, 'Position',[.125 .74 .7 .2])  

 plot(slon,taux_avg,'k')
 hold on
 %plot(slon, taux_trend*12*100, 'g')

 
 cmax = 0.2;
 yrange = (linspace(-1*cmax,cmax,length(redblue)))';
 for j = 1:length(taux_trend)
   
    I_1 = find(yrange <= 12*100*taux_trend(j),1,'last');
    I_2 = find(yrange >= 12*100*taux_trend(j),1,'first');
    
    diff_val = 12*100*taux_trend(j) - yrange(I_1);
    diff_val2 = yrange(I_2)- yrange(I_1);
    frac_val = (diff_val/diff_val2);
    
    taux_colorval = ((1-frac_val)*redblue(I_1,:))+((frac_val)*redblue(I_2,:));
 
    plot(slon(j), 12*100*taux_trend(j),'o','Color',taux_colorval,'MarkerFaceColor',taux_colorval,'MarkerSize',2)
    
 end
 
 plot([slon(1) slon(end)],[0 0],':k')
 
 % Color axis

 caxis([-1*cmax cmax])
colorbar
 

hcb = colorbar;

set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
    'fontweight', 'normal','Position',[0.87 0.748 0.02 0.13])


set(get(hcb,'title'),'string',{'linear trend', '(century^{-1})'},...
    'fontname','arial','fontsize',12,'Position',  [2.5 1.3*cmax 1])

 % Y axis
 
 ylim([-.7 .1])
 
 set(gca,'TickDir','out','Ytick',[-.7 -.5 -.3 -.1 .1])

ylabh = ylabel ('Pressure (N m^{-2})');
  set(ylabh, 'Position', [140 -.3 1])
 
 % X axis
 xlim([150 270])
 set(gca,'TickDir','out','Xtick',[])

 % Border 
set(gca,'box','off')
isholdonque = ishold; 
hold on
ax = axis;
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
plot(ax(1:2),ax(4)*[1,1],'k','linewidth',1)
if isholdonque == 0
hold off
end 

% plot heading
text(150, 0.15,' \tau^x (N m^{-2})','FontSize',12,'fontweight', 'bold');


% island titles
text(168, 0.15,'Maiana','Fontsize',10);
text(180, 0.15,'Howland','Fontsize',10);
text(196, 0.15,'Jarvis','Fontsize',10);

 % islands

plot(200*[1,1],ax(3:4),'color',[0 .6 0]) % Jarvis
plot(184.5*[1,1],ax(3:4),'color',[0 .6 0])  % Howland
plot(173*[1,1],ax(3:4),'color',[0 .6 0]) % Maiana



%

legend ('mean state', 'trend', 'Location', 'SouthEast')
%
%%
hsub = subplot(3,1,2);
set(hsub, 'Position',[.125 .44 .7 .2])
 
 plot(slon, 10^5*we_avg,'k')
 hold on
 %plot(slon, 10^5*we_trend*12*100, 'go')
 plot([slon(1) slon(end)],[0 0],':k')
 

 cmax = 0.5;
 yrange = (linspace(-1*cmax,cmax,length(redblue)))';
 for j = 1:length(10^5*we_trend*12*100)
   
    I_1 = find(yrange <= 10^5*we_trend(j)*12*100,1,'last');
    I_2 = find(yrange >= 10^5*we_trend(j)*12*100,1,'first');
    
    diff_val = 10^5*we_trend(j)*12*100 - yrange(I_1);
    diff_val2 = yrange(I_2)- yrange(I_1);
    frac_val = (diff_val/diff_val2);
    
    we_colorval = ((1-frac_val)*redblue(I_1,:))+((frac_val)*redblue(I_2,:));
 
    plot(slon(j), 10^5*we_trend(j)*12*100,'o','Color',we_colorval,'MarkerFaceColor',we_colorval,'MarkerSize',2)
    
 end

 
 
 
% Color axis

 caxis([-1*cmax cmax])

hcb = colorbar;

set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
    'fontweight', 'normal','Position',[0.87 0.448 0.02 0.13])


set(get(hcb,'title'),'string',{'linear trend', '(century^{-1})'},...
    'fontname','arial','fontsize',12,'Position',  [2.5 1.3*cmax 1])
 
% Y axis
ylim([-1 3])

ylabh = ylabel ('Velocity (x10^{-5} m s^{-1})');
  set(ylabh, 'Position', [140 1 1])


% X axis
xlim([150 270])
set(gca,'TickDir','out','Xtick',[])

 % Border 
set(gca,'box','off')
isholdonque = ishold; 
hold on
ax = axis;
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
plot(ax(1:2),3*[1,1],'k','linewidth',1)
if isholdonque == 0
hold off
end 

% plot heading
text(150, 3.2,'Entrainment Velocity (x10^{-5} m s^{-1})','FontSize',12,'fontweight', 'bold');
% 

% islands

plot(200*[1,1],ax(3:4),'color',[0 .6 0]) % Jarvis
plot(184.5*[1,1],ax(3:4),'color',[0 .6 0])  % Howland
plot(173*[1,1],ax(3:4),'color',[0 .6 0]) % Maiana

 %
hsub =  subplot(3,1,3);
set(hsub, 'Position',[.125 .14 .7 .2])

[a b] = contourf(slon,depth,12*100*u_trend',500);
set(b,'EdgeColor','none')
set(gca,'YDir','reverse')

hold on

[q,h] = contour(slon,depth,u_avg',[-.6,-.3],':k','LineWidth',1);
[qq, hh] = contour(slon,depth,u_avg',[0, 0, 0],'k', 'LineWidth',2);
[qqq, hhh] = contour(slon,depth,u_avg',[0.3, 0.6],'k', 'LineWidth',1);

clabel(q,h, 'manual','fontsize', 10, 'fontweight', 'bold','rotation',0)
clabel(qq,hh, 'manual','fontsize', 10,'fontweight', 'bold','rotation',0)
clabel(qqq,hhh, 'manual','fontsize', 10,'fontweight', 'bold','rotation',0)



 cmax = .3;
 cmin = -.3;
 caxis([cmin cmax])
 
hcb = colorbar;

set(hcb,'YTick',[-1*cmax, 0, cmax],'fontsize', 10,...
    'fontweight', 'normal','Position',[0.87 0.148 0.02 0.13])


set(get(hcb,'title'),'string',{'linear trend', '(century^{-1})'},...
    'fontname','arial','fontsize',12,'Position',  [2.5 1.3*cmax 1])


ylim([4 400]) 

set(gca,'TickDir','out','Ytick',[100 200 300 400])

ylabh = ylabel ('Depth (m)');
  set(ylabh, 'Position', [140 200 1])

 % X axis
xlim([150 270])

set(gca,'TickDir','out','Xtick',[150 180 210 240 270],...
   'XTickLabel',{'150 E', '180', '150 W', '120 W', '90 W'})

 xlabh = xlabel ('Longitude');
 set(xlabh, 'Position', [210 500 1])
 
% Border 
set(gca,'box','off')
isholdonque = ishold; 
hold on
ax = axis;
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
plot(ax(1:2),5*[1,1],'k','linewidth',1)
if isholdonque == 0
hold off
end 

% plot heading
text(150, -15,'Zonal velocity (m s^{-1})','FontSize',12,'fontweight', 'bold');

% islands

plot(200*[1,1],ax(3:4),'color',[0 .6 0]) % Jarvis
plot(184.5*[1,1],ax(3:4),'color',[0 .6 0])  % Howland
plot(173*[1,1],ax(3:4),'color',[0 .6 0]) % Maiana
 