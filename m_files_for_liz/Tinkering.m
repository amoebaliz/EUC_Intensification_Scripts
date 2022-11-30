close all
clear all
clc

lat = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
lon = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
time = ncread('SODA_2.2.6_Trop_Pac_u.cdf','TIME');
depth = ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');
u = ncread('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN');
taux = ncread('SODA_2.2.6_Trop_Pac_taux.cdf','TAUX_ENS_MN');
%%
 % eliminate missing value fill value (-1e+34)
u(u<-10000) = NaN;

 % indexing spatial bounds
 I_lat = find(lat<=2 & lat>=-2);      % For EUC
 I_lon = find(lon>=219 &lon<=221);
 I_dep = find(depth>=10 & depth<=300);
 time0 = 61;
 
 %actual values for EUC
 slat = lat(I_lat);
 slon = lon(I_lon);
 depth = depth(I_dep);

 nt = length(time(time0:end));
 nlon = length(slon);
 nlat = length(slat);
 
 %% Putting time in datenum form

time_2 = zeros(1,length(time(time0:end)));
 
 for j = time0:length(time)
    time_2(j-time0+1) = addtodate(datenum(1865,12,15),double(time(j)),'month');
 end

 %%
u_sub = squeeze(mean(mean(u(I_lon,I_lat,I_dep,time0:end),2),1));
taux_sub = squeeze(mean(mean(taux(I_lon,I_lat,time0:end),2),1));
 
%%
max_filter = 30;
regress_store = ones(max_filter,13)*NaN;
R_store = ones(max_filter,13)*NaN;
signif_store = ones(max_filter,13)*NaN;

for j = 1:max_filter % filtering schemes - error message for j = 18, 21, 25 

   F = filtrage2(max(u_sub)','low', 9, j*12);
   F2 = filtrage2(taux_sub,'low', 9, j*12);

   for k = 1:13 % for whole year and 3 month segments
       if k == 13 % calculate whole year regression
           
           p = polyfit(F2,F,1);
           
           [R P] = corrcoef(F,F2);
         
           regress_store(j,k)=p(1); % store regression slope
           R_store(j,k) = R(1,2);   % store R value
           signif_store(j,k) = P(1,2);  % store P value
       else
           
           if k == 1 % for january centered, leave out 1st year
               i_month = repmat(k-1:k+1,[(nt/12-1) 1]) + repmat(12*(1:(nt/12-1))',[1 3]);
           
           elseif k < 12 && k > 1        
               i_month = repmat(k-1:k+1,[(nt/12) 1]) + repmat(12*(0:(nt/12-1))',[1 3]);
               
           else % for December centered, leave out last year
               i_month = repmat(k-1:k+1,[(nt/12-1) 1]) + repmat(12*(0:(nt/12-2))',[1 3]);
               
           end
           
           F3 = squeeze(mean(F(i_month),2)); % subset F by months
           F4 = squeeze(mean(F2(i_month),2)); % subset F2 by months

           % calculate regresion and calculate significance
           p = polyfit(F4,F3,1);           
           [R P] = corrcoef(F3,F4);        
               
           % store slope and significance value
           
           regress_store(j,k) = p(1);
           R_store(j,k) = R(1,2);   % store R value
           signif_store(j,k) = P(1,2);  % store P value
           plot(F4,F3,'o',F4,p(1)*F4+p(2),'k-')
           %xlim([-.9 -.1])
           %ylim([0.1 0.8])
           pause
           
       end
   end
   
   
   
end

% plot1 = plot(max_taux, max_u,'o');
%%
% for j = 1:nt
%     max_u(j) = max(max(u_sub(:,:,j)));
%     max_taux(j) = min(min(taux_sub(:,j)));
%     %[max_depth(j),~] = find(u_sub(:,j)==max(squeeze(u_sub(:,j))));
%      %sdep(j) = depth(max_depth(j));
%  end
%  

months_store = regress_store(:,1:12);
months_sig = signif_store(:,1:12);

months_store(months_sig>0.05) = NaN;

[c h] = contourf((1:12), (1:max_filter)', months_store,100);
set(h,'edgecolor','none')
load redblue
colormap(redblue)

set(gca,'XTick',[1 3 5 7 9 11],'XTickLabel',['Jan';'Mar';'May';'Jul';'Sep';'Nov'])
caxis([-7 7])

%% 

subplot(3,1,1)

plot(1:length(taux_sub),taux_sub)
title ({'Lon 220'; 'taux'})

subplot(3,1,2)
plot(1:length(taux_sub),max(u_sub))
title ('Max u at lon')

subplot(3,1,3)
plot(1:length(taux_sub),max(u_sub)'./taux_sub)
title('Max u/ taux')

max_euc_sub = max(u_sub);

save taux_EUC_comp taux_sub max_euc_sub
 
%  for m = 0:20
%     
%  if m>0
%      delete(plot1)
%  end
% u_sub = squeeze(mean(mean(u(I_lon+10*m,I_lat,I_dep,time0:end),2),1));
% taux_sub = squeeze(mean(mean(taux(I_lon+10*m,I_lat,time0:end),2),1));
%  
% %% find max zonal velocity 
% max_u = zeros(1,nt);
% max_depth = zeros(1,nt);
% sdep = zeros(1,nt);
% 
% for j = 1:nt
% 
%     [max_depth(j),~,max_u(j)] = find(u_sub(:,j)==max(squeeze(u_sub(:,j))));
%     sdep(j) = depth(max_depth(j));
% end
% 
% plot1 = plot(taux_sub, sdep,'o');
% hold on
% 
% p = polyfit(taux_sub,sdep',1);
% yfit = polyval(p,taux_sub);
% yresid = sdep' - yfit;
% SSresid = sum(yresid.^2);
% SStotal = (length(sdep')-1) * var(sdep');
% rsq = 1 - SSresid/SStotal;
% 
% plot2 = plot(taux_sub,p(1)*taux_sub+p(2),'k-');
% 
% set(gca,'YDir','reverse')
% title(rsq)
% 
% xlim([-2 2])
% ylim([0 300])
% 
% pause
% 
% 
%  end
% 
% % F = filtrage2(max_u,'low',9,7*12);
% % F2 = filtrage2(max_u,'low',9,10*12);
% 
% %save('max_EUC_runmean.mat', 'avgmax_u')
      
%% Calculate the trend

% [u_trend,u_plusminus,u_sig]=trend_stat(max_u,99);
% 
% time_3 = datevec(time_2);
% 
% save ('max.mat', 'time_2', 'u_trend', 'F' , 'max_u')
% 
% figure
% subplot(2,1,1)
% h = plot(time_2,max_u);
% hold on
% z = plot(time_2,F,'k','LineWidth',2);
% datetick('x', 10,'keepticks')
% xlabel('Year')
% ylabel('Meters Per Second (m/s)')
% title('PER = 7*12')
% 
% subplot(2,1,2)
% h = plot(time_2,max_u);
% hold on
% z = plot(time_2,F2,'k','LineWidth',2);
% datetick('x', 10,'keepticks')
% xlabel('Year')
% ylabel('Meters Per Second (m/s)')
% title('PER = 10*12')