%% EUC work
close all

data_sets = {'adcp0n147e_dy.cdf','adcp0n165e_dy.cdf','adcp0n170w_dy.cdf','adcp0n140w_dy.cdf','adcp0n110w_dy.cdf'};
%depth_lev = [4 4 7 10];
%colors = {'b', 'k', 'r', 'b'};

for k = 3:3
    ncid = netcdf.open (char(data_sets(k)),'NC_NOWRITE');
    
    % get variable ID
   
    varid_var = netcdf.inqVarID(ncid,'u_1205');
    depthid = netcdf.inqVarID(ncid,'depth');
    timeid = netcdf.inqVarID(ncid,'time');
   
    U_speed = netcdf.getVar(ncid,varid_var);
    depth = netcdf.getVar(ncid,depthid);
    time = netcdf.getVar(ncid,timeid);
    
    U_speed = squeeze(squeeze(U_speed));
    
    % The "missing value" fill number is 1e+35... so going to remove that

    U_speed_2 =  U_speed;
    U_speed_2(U_speed_2>800) = NaN;

    % get max EUC intensity for every time point and the depth index at
    % which the max occurs
    
    [Y, I] = nanmax(U_speed_2);
    
    EUC_depth = depth(I);
   
    % running mean of both eastward velocitiy and depth
    % choose number of days in mean window:
    n_win = 7*365;
    
    U = runmean(Y',n_win);
    D = runmean(EUC_depth,n_win);
    
    [gtime] = gregorian(double(time));
    x_date = datenum(gtime);
    
    % plotting... same thing we do every night, Pinky
    
    a = datenum(1988,1,1);
    b = datenum(2012,1,1);
    
    subplot(2,5,k)
    plot(x_date,U,'LineWidth',1)
    xlim([a b])
    datetick('x','yyyy','keeplimits')
    
    subplot(2,5,k+5)
    plot(x_date,D,'LineWidth',1)
    axis([a b 10 200])
    set(gca,'YDir','reverse')
    datetick('x','yyyy','keeplimits')
    
    % plot(x_date,T,char(colors(k)),'LineWidth',1)
    %plot(x_date,T/abs(max(T)),char(colors(k)),'LineWidth',1)
    %hold on
end


