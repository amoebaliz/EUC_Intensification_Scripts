close all 
clear all
clc
 ncid = netcdf.open ('u_eqpac_SODA_extended.nc','NC_NOWRITE'); 
    
    % get variable IDs for SEC
    varid_lat = netcdf.inqVarID(ncid,'lat');
    varid_lon = netcdf.inqVarID(ncid,'lon');
    varid_time = netcdf.inqVarID(ncid,'time');
    varid_depth = netcdf.inqVarID(ncid,'depth');
    
    lat = netcdf.getVar(ncid,varid_lat);
    lon = netcdf.getVar(ncid,varid_lon);
    depth = netcdf.getVar(ncid,varid_depth);
    time = netcdf.getVar(ncid,varid_time);
    
    % get variable IDs for EUC
    varid_u_euc = netcdf.inqVarID(ncid,'u');
    u_euc = netcdf.getVar(ncid,varid_u_euc);
    
    %% Indexing spatial bounds
 
 I_lat_1 = find(lat<=2 & lat>=-2);      % Gilberts
 I_lat_2 = find(lat<=1 & lat>=0);      % Howland & Baker
 I_lat_3 = find(lat<=0 & lat>=-0.5);      % Jarvis
 
 
 I_lon_1 = find(lon>172 &lon<175);      % Gilberts
 I_lon_2 = find(lon>183 &lon<184);      % Howland & Baker
 I_lon_3 = find(lon>200 &lon<=201);      % Jarvis
 
 I_depth = find(depth>=0 &depth<=500); 
 
 %% 
 
 EUC_1 = u_euc(I_lon_1,I_lat_1, I_depth,:);
 EUC_2 = u_euc(I_lon_2,I_lat_2, I_depth,:);
 EUC_3 = u_euc(I_lon_3,I_lat_3, I_depth,:);
 
%% grid "height" vector
height = zeros(length(I_depth),1);  % EUC
depth_int = diff([0;depth(I_depth)]);

for j = 2:length(I_depth)
    height(j-1) = 0.5*(depth_int(j) + depth_int(j-1));  
end
height(1) = height(1)+ 0.5*depth_int(1);
height(end) = depth_int(end);

%% width of grid box
y_dim = 0.5*110.574*100;
win = 13;

euc_trend = 50*ones(1,3); euc_plusminus = euc_trend; euc_sig = euc_trend;
mean_euc_1D = zeros(1656,3);
t= 0;
for j = 1:3

% EUC

I_lat = eval(['I_lat_',num2str(j)]);

euc_grid_area_1 = repmat(height*y_dim,1,numel(find(I_lat)));
euc_grid_area_2 = repmat(euc_grid_area_1,[1 1 1]);
euc_grid_area_3 = permute(euc_grid_area_2,[3 2 1]);
euc_grid_area_4 = repmat(euc_grid_area_3,[1 1 1 length(time)]);

% multiply by data

u_euc_island = eval(['EUC_',num2str(j)]);

if j ==3
    u_euc_island = (mean(u_euc_island(:,:,:,:),1));
else
    u_euc_island = squeeze(mean(u_euc_island(:,:,:,:),1));
end
u_euc_island(u_euc_island<=0) = 0;
u_euc_island(isnan(u_euc_island)) = 0;

if j == 3
    transport = (u_euc_island.*euc_grid_area_4);
else
transport = (u_euc_island.*squeeze(euc_grid_area_4));
end

if j == 3
    trans_sum = squeeze(sum(transport(:,:,:,:),3));
else
    trans_sum = sum(sum(transport(:,:,:),1),2);
end
u_euc_isnan(u_euc_island<=0) = 0;
area_calc = sum(sum(sum(euc_grid_area_4(:,:,:,1),1),2),3);

euc_1D = trans_sum/area_calc;               % should be left with number over all time
size(euc_1D)
 % get running mean       
 
mean_euc_1D(:,j) = runmean(squeeze(euc_1D),win);
 
 %calculate trend statistics
 [euc_trend(j),euc_plusminus(j),euc_sig(j)]=trend_stat(mean_euc_1D(:,j),99);
 t= t+1
end
 
save('Island_data.mat', 'mean_euc_1D', 'euc_trend', 'euc_plusminus', 'euc_sig')