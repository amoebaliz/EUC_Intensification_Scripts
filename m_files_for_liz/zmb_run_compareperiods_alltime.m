% This is a script to run the isopycnal zonal momentum budget on REDONDA
% and compare quarter-record-length means.

% set up workspace

close all, clear all
%addpath('/data1/Code/MATLAB/')
%addpath('/data1/Code/MATLAB/seawater_ver3_3/')
load redblue.mat
time0 = 61;
 q1 = time0:time0+35*12-1; 
 q4 = 1717-35*12:1716;
 i_mon_mam = zeros(35,3);
 i_mon_jja = zeros(35,3);

 for i_nyr = 1:35;
                    
    i_mon_mam(i_nyr,:) = (3:5)+12*(i_nyr-1);
    i_mon_jja(i_nyr,:) = (6:8)+12*(i_nyr-1);
    
 end
i_mon_mam = sort(i_mon_mam(:));
i_mon_jja = sort(i_mon_jja(:));
                
% specify SODA quarters

% Create mean state (1) or seasonal (2) data sets

% open full SODA record variables

            lon=ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
            lat=ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
            depth=ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');

%% u

u=ncread('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN');

u1=u(:,:,:,q1);      
u4=u(:,:,:,q4);
clear u

u1_mam = u1(:,:,:,i_mon_mam);
u4_mam = u4(:,:,:,i_mon_mam);

u1_jja = u1(:,:,:,i_mon_jja);
u4_jja = u4(:,:,:,i_mon_jja);

%% v

v=ncread('SODA_2.2.6_Trop_Pac_v.cdf','V_ENS_MN');
v1=v(:,:,:,q1);      
v4=v(:,:,:,q4);
clear v

v1_mam = v1(:,:,:,i_mon_mam);
v4_mam = v4(:,:,:,i_mon_mam);

v1_jja = v1(:,:,:,i_mon_jja);
v4_jja = v4(:,:,:,i_mon_jja);
                
%% w

w=ncread('SODA_2.2.6_Trop_Pac_w.cdf','W_ENS_MN');
w1=w(:,:,:,q1);      
w4=w(:,:,:,q4);
clear w

w1_mam = w1(:,:,:,i_mon_mam);
w4_mam = w4(:,:,:,i_mon_mam);

w1_jja = w1(:,:,:,i_mon_jja);
w4_jja = w4(:,:,:,i_mon_jja);

%% temp

temp=ncread('SODA_2.2.6_Trop_Pac_temp.cdf','TEMP_ENS_MN');
temp1=temp(:,:,:,q1);      
temp4=temp(:,:,:,q4);
clear temp

temp1_mam = temp1(:,:,:,i_mon_mam);
temp4_mam = temp4(:,:,:,i_mon_mam);

temp1_jja = temp1(:,:,:,i_mon_jja);
temp4_jja = temp4(:,:,:,i_mon_jja);

%% salt

salt=ncread('SODA_2.2.6_Trop_Pac_salt.cdf','SALT_ENS_MN');
salt1=salt(:,:,:,q1);      
salt4=salt(:,:,:,q4);
clear salt

salt1_mam = salt1(:,:,:,i_mon_mam);
salt4_mam = salt4(:,:,:,i_mon_mam);

salt1_jja = salt1(:,:,:,i_mon_jja);
salt4_jja = salt4(:,:,:,i_mon_jja);

%% ssh

ssh=ncread('SODA_2.2.6_Trop_Pac_ssh.cdf','SSH_ENS_MN');
ssh1=ssh(:,:,q1);      
ssh4=ssh(:,:,q4);
clear ssh

ssh1_mam = ssh1(:,:,i_mon_mam);
ssh4_mam = ssh4(:,:,i_mon_mam);

ssh1_jja = ssh1(:,:,i_mon_jja);
ssh4_jja = ssh4(:,:,i_mon_jja);
 
%% taux
taux=ncread('SODA_2.2.6_Trop_Pac_taux.cdf','TAUX_ENS_MN');
taux1=taux(:,:,q1);      
taux4=taux(:,:,q4);
clear taux

taux1_mam = taux1(:,:,i_mon_mam);
taux4_mam = taux4(:,:,i_mon_mam);

taux1_jja = taux1(:,:,i_mon_jja);
taux4_jja = taux4(:,:,i_mon_jja);

%%% run ZMB for each quarter-record-length seasonal mean %%%


%% q1 mam

for t = 1:35*3
    
    [~,latx,depth,~,~,ududx,vdudy,wdudz,zpgf,~,frich,fricv] = ...
        zmb_isopycnal(lon,lat,depth,squeeze(u1_mam(:,:,:,t)),squeeze(v1_mam(:,:,:,t)),...
        squeeze(w1_mam(:,:,:,t)),squeeze(temp1_mam(:,:,:,t)),squeeze(salt1_mam(:,:,:,t)),...
        squeeze(ssh1_mam(:,:,t)));
    
    ududx_store(:,:,:,t) = ududx;
    vdudy_store(:,:,:,t) = vdudy;
    wdudx_store(:,:,:,t) = wdudz;
    
    zpgf_store(:,:,:,t) = zpgf;
    frich_store(:,:,:,t) = frich;
    fricv_store(:,:,:,t) = fricv;
   
    disp(t)
end

save MAM_quarter1_sens lon latx depth taux1_mam ssh1_mam ududx_store vdudy_store wdudz_store zpgf_store frich_store fricv_store

%% q4 mam
% 
% for t = 1:35*3
%     
%     [~,latx(:,t),depth(:,t),~,~,ududx(:,:,:,t),vdudy(:,:,:,t),wdudz(:,:,:,t),zpgf(:,:,:,t),~,frich(:,:,:,t),fricv(:,:,:,t)] = ...
%         zmb_isopycnal(lon,lat,depth,squeeze(u1_mam(:,:,:,t)),squeeze(v1_mam(:,:,:,t)),...
%         squeeze(w1_mam(:,:,:,t)),squeeze(temp1_mam(:,:,t)),squeeze(salt1_mam(:,:,t)),...
%         squeeze(ssh1_mam(:,:,t)));
% end
% 
% 

