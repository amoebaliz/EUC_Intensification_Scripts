% EUC SODA data

close all
    ncid = netcdf.open ('u_eqpac_SODA.nc','NC_NOWRITE');
    
     % lon x lat x depth x time

    % get variable ID
    varid_lat = netcdf.inqVarID(ncid,'lat');
    varid_lon = netcdf.inqVarID(ncid,'lon');
    varid_time = netcdf.inqVarID(ncid,'time');
    varid_u = netcdf.inqVarID(ncid,'u');
    varid_depth = netcdf.inqVarID(ncid,'depth');
    
    % get variables
    time = netcdf.getVar(ncid,varid_time);
    lat = netcdf.getVar(ncid,varid_lat);
    lon = netcdf.getVar(ncid,varid_lon);
    u = netcdf.getVar(ncid,varid_u);
    depth = netcdf.getVar(ncid,varid_depth);
    
 % eliminate missing value fill value
 u(u<-10) = NaN;
 
  % putting time in datenum format
 time_2 = datenum(1960,(1+double(time)),1);
 
 I_lat = find(lat<1 & lat>-1);          % On Equator
 % I_lon = find(lon>189.5 & lon<190.5);   % At comparable latitude
 % I_lon = find(lon>219.5 & lon<220.5);
 I_lon = find(lon>249.5 & lon<250.5);
 
 %I_dep = find(depth>40);si
 
 max_u = squeeze(squeeze(squeeze(max(mean(mean(u(I_lon,I_lat,:,:)))))));
 
 %% getting ADCP data
 
    % Files = {'adcp0n110w_mon.cdf','adcp0n140w_mon.cdf','adcp0n170w_mon.cdf'}
    % 110w time = days since 1991-05-17
    % 140w time = days since 1990-05-01
    % 170w time = days since 1988-05-16
    
    ADCPncid = netcdf.open ('adcp0n110w_mon.cdf','NC_NOWRITE');
 
    ADCPid_lat = netcdf.inqVarID(ADCPncid,'lat');
    ADCPid_lon = netcdf.inqVarID(ADCPncid,'lon');
    ADCPid_time = netcdf.inqVarID(ADCPncid,'time');
    ADCPid_u = netcdf.inqVarID(ADCPncid,'u_1205');
    ADCPid_depth = netcdf.inqVarID(ADCPncid,'depth');
 
    % get variables
    ADCPtime = netcdf.getVar(ADCPncid,ADCPid_time);
    ADCPlat = netcdf.getVar(ADCPncid,ADCPid_lat);
    ADCPlon = netcdf.getVar(ADCPncid,ADCPid_lon);
    ADCPu = netcdf.getVar(ADCPncid,ADCPid_u);
    ADCPdepth = netcdf.getVar(ADCPncid,ADCPid_depth);
    
    % eliminate missing value fill value
    ADCPu(ADCPu>500) = NaN;
    
    %convert to m/s
    ADCPu = ADCPu/100;
    
    % putting time in datenum format
    time_3 = datenum(1991,5,(17+double(ADCPtime)));
    
    max_ADCPu = squeeze(squeeze(max(ADCPu)));
 %% running means
 
 win = 12+1;    %   
 avgmax_u = runmean(max_u,win);
 
 win2 = 12+1;    %   
 avgmax_ADCPu = runmean(max_ADCPu,win2);
 
 %%
 
%  slon = lon(I_lon);
%  slat = lat(I_lat);
%  sdep = depth(I_dep);
%  
%  nt = length(time);
%  nlon = length(slon);
%  nlat = length(slat);
%  ndep = length(sdep);
%  
%  Istore = [];
%  Jstore = [];
%  Kstore = [];
%  count = [];
%  
%  for jj = 1:nt
%     
%     
%     [I J K] = ind2sub([nlon nlat ndep],find(u(I_lon,I_lat,I_dep,jj) == squeeze(squeeze(squeeze(max(max(max(u(I_lon,I_lat,I_dep,jj)))))))));
%     %[I J K] = ind2sub([nlon nlat ndep],find(u(:,:,:,jj) == squeeze(squeeze(squeeze(max(max(max(u(:,:,:,jj)))))))));
%     
% %     I = find(u == squeeze(squeeze(squeeze(max(max(max(u(:,:,:,jj))))))));
%     count(jj)=numel(I);
%     
%     Istore = [Istore; I];
%     Jstore = [Jstore; J];
%     Kstore = [Kstore; K];
%  end
%  
%  close all
%  plot3(slon(Istore),slat(Jstore),sdep(Kstore),'o')
%   set(gca,'ZDir','reverse')
 
 %% plotting the data
 
 plot(time_2,avgmax_u,time_3,avgmax_ADCPu)

 datetick('x','yyyy','keeplimits')

 legend('SODA', 'ADCP', 'Location', 'northwest')