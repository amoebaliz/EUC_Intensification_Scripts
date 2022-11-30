% ADCP comparison

clear all
close all
clc

lon_soda = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
lat_soda = ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
time_soda = ncread('SODA_2.2.6_Trop_Pac_u.cdf','TIME');
depth_soda = ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');
u_soda = ncread('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN');
u_soda(u_soda<-1000)=NaN;

%%%% Time info for soda
%
time_soda_2 = zeros(1,length(time_soda));
 for j = 1:length(time_soda)
    time_soda_2(j) = addtodate(datenum(1865,12,15),double(time_soda(j)),'month');
 end

% Names of TAO buoy data
tao_sites = {'adcp0n165e_dy.cdf','adcp0n170w_dy.cdf',...
    'adcp0n140w_dy.cdf', 'adcp0n110w_dy.cdf'};
%%
store_R = zeros(20,4);
stored_data = cell(4,4);
%%
for j = 1:4
    
    
    time_att = ncreadatt(char(tao_sites(j)),'time','units');
    time = ncread(char(tao_sites(j)),'time');
    depth = ncread(char(tao_sites(j)),'depth');
    lon = ncread(char(tao_sites(j)),'lon');
    u = ncread(char(tao_sites(j)),'u_1205');
    
    u = squeeze(u/100); %convert to m/s
    u(u>1000)=NaN;
    
    %%%%%%%%%  Regrid ADCP data to match the SODA structure  %%%%%%%%%%
    u_adcp_regrid = zeros(length(depth_soda),length(time));
    warning('off','all');
    for z= 1:length(time)
                
        u_adcp_regrid(:,z)=interp1(depth,squeeze(u(:,z)),depth_soda);

    end; clear z
    warning('on','all');    
    
    
%%
     %%%%%%%%  Get Monthly Averages for TAO ADCP Data   %%%%%%%%%%%
    % 
    % Putting time in datenum form
     time_2 = zeros(1,length(time));  
     timen = time;
        for k = 1:length(time)
            time_2(k) = addtodate(datenum(str2num(time_att(12:15)),...
            str2num(time_att(17:18)),str2num(time_att(20:21))),...
            double(timen(k)),'day');
        end
   
    % setting up to get monthly average using grpstats
     date_nums = datevec(time_2);   
     T = [date_nums(:,1:2) zeros(length(date_nums),1)];
     
     
     %initialize avg depth vriable
     % U_avg_depth = zeros(length(depth_soda),length(T_2));
    clear U_avg_depth
    
    for k = 1:length(depth_soda)
        T(:,3) = u_adcp_regrid(k,:);
        T_2 = grpstats(T,T(:,1:2));
        U_avg_depth(k,:)= T_2(:,3)';
    end

      % getting dates back
      ADCP_date_nums = zeros(length(T_2),1);
      
    for kk = 1:length(T_2)
        ADCP_date_nums(kk) = datenum(T_2(kk,1),T_2(kk,2),15);
    end
    
    ADCP_date_nums= ADCP_date_nums';
    
    % Identify where the corresponding longitude occurs in SODA
   
    i_lon = find(lon_soda< lon+.5 & lon_soda> lon-.5);
    %%
    % find max value through depth
    
    test_nan = mean(mean(u_soda(i_lon,10:11,:,:),1),2);
    
    [max_soda_series, i_soda_max]= max(mean(mean(u_soda(i_lon,10:11,:,:),1),2));
    max_soda_series = squeeze(max_soda_series);
    i_soda_max = squeeze(i_soda_max);
    
    % find where in SODA does the time comparison start
    i_soda_start = find(time_soda_2 == ADCP_date_nums(1));
    
    [max_u_adcp,i_adcp_max] = max(U_avg_depth);
    max_u_adcp = squeeze(max_u_adcp);
    
    max_adcp_depths = depth_soda(i_adcp_max);
    max_soda_depths = depth_soda(i_soda_max);
    
    max_soda_series(max_soda_depths>300) = NaN;
    max_soda_depths(max_soda_depths>300) = NaN;
%%    

if length(max_u_adcp)>length(max_soda_series(i_soda_start:end))
        
        u_adcp = U_avg_depth(:,1:length(max_soda_series(i_soda_start:end)));
        u_soda_2 = squeeze(mean(mean(u_soda(i_lon,10:11,:,i_soda_start:end),2),1));
        
else
        u_adcp = U_avg_depth(:,1:end);
        u_soda_2 = squeeze(mean(mean(u_soda(i_lon,10:11,:,i_soda_start:(i_soda_start+length(max_u_adcp)-1)),2),1));
        
end

velocity_R = corrcoef(u_adcp(:),u_soda_2(:),'rows','pairwise');
%%

for teal = 1:20
velocity_R_2 = corrcoef(u_adcp(teal,:),u_soda_2(teal,:),'rows','pairwise');
store_R(teal,j) = velocity_R_2(1,2);
end
%% Difference in plots: 
figure(3)
     subplot(1,length(tao_sites),j)
     %plot ADCP
      m = plot(1:length(u_adcp(:)),u_adcp(:),'k');

      hold on
      %plot SODA
      mm = plot(1:length(u_soda_2(:)),u_soda_2(:),'b-.');        
    
      legend('ADCP','SODA')         
         
        if j == 1
            ylabel({'Velocity','meters/second'})

        end
        
        title(lon)
        ylim([-2 2])
        xlim([1 length(u_soda_2(:))])
        text(1500,-1.5,['R=' num2str(velocity_R(1,2))])
            

    if length(max_u_adcp)>length(max_soda_series(i_soda_start:end))
        velocity_R = corrcoef(max_u_adcp(1:length(max_soda_series(i_soda_start:end)))',max_soda_series(i_soda_start:end),'rows','pairwise');
        depth_R = corrcoef(max_adcp_depths(1:length(max_soda_series(i_soda_start:end)))',max_soda_depths(i_soda_start:end),'rows','pairwise');
        bias_vals = nanmean(max_adcp_depths(1:length(max_soda_series(i_soda_start:end)))-max_soda_depths(i_soda_start:end));
    else
        velocity_R = corrcoef(max_u_adcp', max_soda_series(i_soda_start:i_soda_start+length(max_u_adcp)-1),'rows','pairwise');
        depth_R = corrcoef(max_adcp_depths,max_soda_depths(i_soda_start:i_soda_start+length(max_u_adcp)-1),'rows','pairwise');
        bias_vals = nanmean(max_adcp_depths - max_soda_depths(i_soda_start:i_soda_start+length(max_u_adcp)-1));
    end
%%    
    figure(4)
     subplot(1,length(tao_sites),j)
        %plot ADCP
        plot(ADCP_date_nums,max_u_adcp,'k')
        datetick
        hold on
        %plot SODA
        plot(time_soda_2(i_soda_start:end),max_soda_series(i_soda_start:end),'b-.')
        
        if j == 1
            ylabel({'Max Velocity','meters/second'})
        end
       
        title(lon)
        ylim([0 2])
       if j<3
            text(time_soda_2(floor(i_soda_start + (length(time_soda_2(i_soda_start:end))/1.5))),1.6,['R=', num2str(velocity_R(1,2))])
       else
            text(time_soda_2(floor(i_soda_start + (length(time_soda_2(i_soda_start:end))/1.5))),0.2,['R=', num2str(velocity_R(1,2))])
       end
            legend('ADCP','SODA')
 %%      
     stored_data{1,j} = datevec(ADCP_date_nums);
     stored_data{2,j} = max_u_adcp;
     stored_data{3,j} = max_soda_series(i_soda_start:end);
     stored_data{4,j} = velocity_R(1,2);
            
            
            
            
            
            
            
            
        subplot(2,length(tao_sites),j+4)
        plot(ADCP_date_nums,max_adcp_depths,'k')
        datetick
        hold on
        plot(time_soda_2(i_soda_start:end),max_soda_depths(i_soda_start:end),'b-.')
        
        if j == 1
            ylabel({'Depth of Max Velocity','meters'})
        end
        legend('ADCP','SODA')
       
        text(time_soda_2(i_soda_start+(length(time_soda_2(i_soda_start:end))/2)),400,['R=' num2str(depth_R(1,2))])

        ylim([0 500])
        set(gca,'YDir','reverse')
        
%%% Figure comparing the time series on SODA and ADCP    
%     
%     figure(1)
% set(gcf, 'Renderer', 'painter')
%     
% subplot(4,2,2*(j-1)+1)
% m = pcolor(repmat(time_soda_2(i_soda_start:end),[length(depth_soda),1])',repmat(depth_soda',...
%     [length(time_soda_2(i_soda_start:end)),1]),...
%     double(squeeze(mean(mean(u_soda(i_lon,10:11,:,i_soda_start:end),1),2))'));
% 
% xlim([datenum(1992,1,1) datenum(2008,1,1)])   
% set(gca,'YDir','reverse')
% set(m,'EdgeColor','none');
% ylabel('Depth (m)')
% caxis([-1.5 1.5])
% ylim([0 300])    
% datetick('x','keeplimits')
% 
% subplot(4,2,2*j)
% m = pcolor(repmat(ADCP_date_nums',[length(depth_soda),1])',... 
%     repmat(depth_soda',[length(ADCP_date_nums),1]),double( U_avg_depth'));
% set(gca,'YDir','reverse')
% set(m,'EdgeColor','none');
% xlim([datenum(1992,1,1) datenum(2008,1,1)])   
% 
% ylabel('Depth (m)')
% caxis([-1.5 1.5])
% ylim([0 300])
% datetick('x','keeplimits')
    
% 
%     
%     clear('ADCP_date_nums')
%     clear('U_avg_depth')

%checking time series - regridded vs. original ADCP
% figure(1)
% subplot(4,2,2*(j-1)+1)
% m = pcolor(repmat(time',[length(depth),1])',repmat(depth',[length(time),1]),double(u'));
% set(gca,'YDir','reverse')
% set(m,'EdgeColor','none');
% ylabel('Depth (m)')
% caxis([-1.5 1.5])
% ylim([0 300])    
%     
% subplot(4,2,2*j)
% m = pcolor(repmat(time',[length(depth_soda),1])',repmat(depth_soda',[length(time),1]),double(u_adcp_regrid'));
% set(gca,'YDir','reverse')
% set(m,'EdgeColor','none');
% 
% ylabel('Depth (m)')
% caxis([-1.5 1.5])
% ylim([0 300])
% 
% %bias_val = mean(max_adcp_depths-max_soda_series(i_soda_start:end));
% disp(bias_vals)
end
% 
%%
% close all
% 
% I_lon = find(lon_soda >= 150 & lon_soda <= 270);
% meep = mean(u_soda(I_lon, 10:11,:,i_soda_start:end),2);
% 
% %nan_index = find(isnan(max_soda_series(i_soda_start:end)));
% nan_index = [45, 48, 74, 82];
% 
% figure
% load redblue
% for j = 1:4
%     
%    k = nan_index(j);
%    
%     h = pcolor(lon_soda(I_lon),depth_soda,squeeze(meep(:,:,:,k))');
%     set(h,'EdgeColor','none')
%     set(gca,'YDir','reverse')
%     
%     hold on
%     plot([mean(lon_soda(i_lon)) mean(lon_soda(i_lon))],[depth_soda(1) depth_soda(end)])
%     colormap(redblue)
%     caxis([-1 1])
%     
%     MON = datevec(ADCP_date_nums(k));
%     
%     title([num2str(MON(2)),' ,',num2str(MON(1))])
%     pause
%     
% end
