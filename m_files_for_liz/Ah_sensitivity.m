% This is a script to run the isopycnal zonal momentum budget on REDONDA
% and compare quarter-record-length means.

% set up workspace

close all, clear all
%addpath('/data1/Code/MATLAB/')
%addpath('/data1/Code/MATLAB/seawater_ver3_3/')
load redblue.mat
time0 = 61;
% specify SODA quarters

% Create mean state (1) or seasonal (2) data sets
for get_mean = 1:3

    % Prescrive the season: 0 = mean state; 1 = MAM; 2 = JJA
   
    if get_mean == 2
    
        seas = 3:5; s=1;      % Mar/Apr/May
    
    elseif get_mean == 3
     
        seas = 6:8; s=2;      % June/July/Aug
        
    end
  
             q1 = time0:time0+35*12; 
             %q2 = 433:864; 
             %q3 = 853:1284 ; 
             q4 = 1717-35*12:1716;
  

            % open full SODA record variables

            lon=ncread('SODA_2.2.6_Trop_Pac_u.cdf','LON241_560');
            lat=ncread('SODA_2.2.6_Trop_Pac_u.cdf','LAT142_161');
            depth=ncread('SODA_2.2.6_Trop_Pac_u.cdf','DEPTH1_20');

            %
            u=ncread('SODA_2.2.6_Trop_Pac_u.cdf','U_ENS_MN');

            if get_mean==1
                u_mean = mean(u(:,:,:,time0:end),4);
            else
                  u=permute(u,[4 3 2 1]);         % permute for use in clim code: Time x depth x lat x lon

                 u1=clim3d(u(q1,:,:,:));      % compute quarter-record-length clims
                 u4=clim3d(u(q4,:,:,:));
                % u2=clim3d(u(q2,:,:,:));      
                % u3=clim3d(u(q3,:,:,:));

                u1=permute(u1,[4 3 2 1]);       % permute back for use in zmb code: lon x lat x depth x time
                u4=permute(u4,[4 3 2 1]);
                % u2=permute(u2,[4 3 2 1]);       
                % u3=permute(u3,[4 3 2 1]);

                u1=mean(u1(:,:,:,seas),4);      % specify/form climped seasonal mean
                u4=mean(u4(:,:,:,seas),4);
                % u2=mean(u2(:,:,:,seas),4);      
                % u3=mean(u3(:,:,:,seas),4);
            end

            clear u

            v=ncread('SODA_2.2.6_Trop_Pac_v.cdf','V_ENS_MN');

            if get_mean == 1
              v_mean = mean(v(:,:,:,time0:end),4);  
            else
                 v=permute(v,[4 3 2 1]);         % permute for use in clim code: Time x depth x lat x lon
                 v1=clim3d(v(q1,:,:,:));      % compute quarter-record-length clims
                 v4=clim3d(v(q4,:,:,:));
                % v2=clim3d(v(q2,:,:,:));      
                % v3=clim3d(v(q3,:,:,:));

                v1=permute(v1,[4 3 2 1]);       % permute back for use in zmb code: lon x lat x depth x time
                v4=permute(v4,[4 3 2 1]);
                % v2=permute(v2,[4 3 2 1]);       
                % v3=permute(v3,[4 3 2 1]);

                v1=mean(v1(:,:,:,seas),4);      % specify/form climped seasonal mean
                v4=mean(v4(:,:,:,seas),4);
                % v2=mean(v2(:,:,:,seas),4);      % specify/form climped seasonal mean
                % v3=mean(v3(:,:,:,seas),4);

            end

            clear v

            w=ncread('SODA_2.2.6_Trop_Pac_w.cdf','W_ENS_MN');

            if get_mean == 1
                w_mean = mean(w(:,:,:,time0:end),4);
            else

            w=permute(w,[4 3 2 1]);         % permute for use in clim code: Time x depth x lat x lon
            w1=clim3d(w(q1,:,:,:));      % compute quarter-record-length clims
            w4=clim3d(w(q4,:,:,:));
            % w2=clim3d(w(q2,:,:,:));      
            % w3=clim3d(w(q3,:,:,:));

            w1=permute(w1,[4 3 2 1]);       % permute back for use in zmb code: lon x lat x depth x time
            w4=permute(w4,[4 3 2 1]);
            % w2=permute(w2,[4 3 2 1]);       
            % w3=permute(w3,[4 3 2 1]);

            w1=mean(w1(:,:,:,seas),4);      % specify/form climped seasonal mean
            w4=mean(w4(:,:,:,seas),4);
            % w2=mean(w2(:,:,:,seas),4);      
            % w3=mean(w3(:,:,:,seas),4);

            end

            clear w

            temp=ncread('SODA_2.2.6_Trop_Pac_temp.cdf','TEMP_ENS_MN');

            if get_mean == 1

                temp_mean = mean(temp(:,:,:,time0:end),4);
            else
                temp=permute(temp,[4 3 2 1]);   % permute for use in clim code: Time x depth x lat x lon
                temp1=clim3d(temp(q1,:,:,:));% compute quarter-record-length clims
                temp4=clim3d(temp(q4,:,:,:));
                % temp2=clim3d(temp(q2,:,:,:));
                % temp3=clim3d(temp(q3,:,:,:));

                temp1=permute(temp1,[4 3 2 1]); % permute back for use in zmb code: lon x lat x depth x time
                temp4=permute(temp4,[4 3 2 1]);
                % temp2=permute(temp2,[4 3 2 1]); 
                % temp3=permute(temp3,[4 3 2 1]);

                temp1=mean(temp1(:,:,:,seas),4);% specify/form climped seasonal mean
                temp4=mean(temp4(:,:,:,seas),4);
                % temp2=mean(temp2(:,:,:,seas),4);
                % temp3=mean(temp3(:,:,:,seas),4);

            end

            clear temp

            salt=ncread('SODA_2.2.6_Trop_Pac_salt.cdf','SALT_ENS_MN');

            if get_mean == 1
                salt_mean = mean(salt(:,:,:,time0:end),4);
            else
                salt=permute(salt,[4 3 2 1]);   % permute for use in clim code: Time x depth x lat x lon
                salt1=clim3d(salt(q1,:,:,:));% compute quarter-record-length clims
                salt4=clim3d(salt(q4,:,:,:));
                % salt2=clim3d(salt(q2,:,:,:));
                % salt3=clim3d(salt(q3,:,:,:));

                salt1=permute(salt1,[4 3 2 1]); % permute back for use in zmb code: lon x lat x depth x time
                salt4=permute(salt4,[4 3 2 1]);
                % salt2=permute(salt2,[4 3 2 1]); 
                % salt3=permute(salt3,[4 3 2 1]);

                salt1=mean(salt1(:,:,:,seas),4);% specify/form climped seasonal mean   
                salt4=mean(salt4(:,:,:,seas),4);
                % salt2=mean(salt2(:,:,:,seas),4);
                % salt3=mean(salt3(:,:,:,seas),4);
            end

            clear salt

            ssh=ncread('SODA_2.2.6_Trop_Pac_ssh.cdf','SSH_ENS_MN');

            if get_mean == 1
                ssh_mean = mean(ssh(:,:,time0:end),3);

            else

                ssh=permute(ssh,[3 2 1]);       % permute for use in clim code: Time x depth x lat x lon
                ssh1=clim2d(ssh(q1,:,:));    % compute quarter-record-length clims
                ssh4=clim2d(ssh(q4,:,:));
                %  ssh2=clim2d(ssh(q2,:,:));    
                %  ssh3=clim2d(ssh(q3,:,:));

                ssh1=permute(ssh1,[3 2 1]);     % permute back for use in zmb code: lon x lat x depth x time
                ssh4=permute(ssh4,[3 2 1]);
                % ssh2=permute(ssh2,[3 2 1]);     
                % ssh3=permute(ssh3,[3 2 1]);

                ssh1=mean(ssh1(:,:,seas),3);    % specify/form climped seasonal mean
                ssh4=mean(ssh4(:,:,seas),3);
                % ssh2=mean(ssh2(:,:,seas),3);    
                % ssh3=mean(ssh3(:,:,seas),3);

            end

            clear ssh

            taux=ncread('SODA_2.2.6_Trop_Pac_taux.cdf','TAUX_ENS_MN');

            if get_mean == 1

                taux_mean = mean(taux(:,:,time0:end),3);

            else

            taux=permute(taux,[3 2 1]); 
            taux1=clim2d(taux(q1,:,:));    % compute quarter-record-length clims
            taux4=clim2d(taux(q4,:,:));
            % taux2=clim2d(taux(q2,:,:));    
            % taux3=clim2d(taux(q3,:,:));

            taux1=permute(taux1,[3 2 1]);     % permute back for use in zmb code: lon x lat x depth x time
            taux4=permute(taux4,[3 2 1]);
            % taux2=permute(taux2,[3 2 1]);     
            % taux3=permute(taux3,[3 2 1]);

            taux1=squeeze(mean(taux1(:,:,seas),3));    % specify/form climped seasonal mean
            taux4=squeeze(mean(taux4(:,:,seas),3));
            % taux2=squeeze(mean(taux2(:,:,seas),3));    
            % taux3=squeeze(mean(taux3(:,:,seas),3));

            end

            clear taux

            %%
            % run ZMB for each quarter-record-length seasonal mean
                    if get_mean == 1
                         [~,latx,depth5,temp5,u5,dudt5,ududx5,vdudy5,wdudz5,zpgf5,corf5,frich5,fricv5] = zmb_isopycnal(lon,lat,depth,u_mean,v_mean,w_mean,temp_mean,salt_mean,ssh_mean);

                            save mean_fields_8000 lon latx depth5 u5 ssh_mean taux_mean... 
                                 dudt5 ududx5 vdudy5 wdudz5 zpgf5 corf5 frich5 fricv5

                     else
                         [~,~,depth1,temp1,u1,dudt1,ududx1,vdudy1,wdudz1,zpgf1,corf1,frich1,fricv1] = zmb_isopycnal(lon,lat,depth,u1,v1,w1,temp1,salt1,ssh1);
                         [~,latx,depth4,temp4,u4,dudt4,ududx4,vdudy4,wdudz4,zpgf4,corf4,frich4,fricv4] = zmb_isopycnal(lon,lat,depth,u4,v4,w4,temp4,salt4,ssh4);
                        % 
                        % [~,latx,depth2,u2,dudt2,ududx2,vdudy2,wdudz2,zpgf2,corf2,frich2,fricv2] = zmb_isopycnal(lon,lat,depth,u2,v2,w2,temp2,salt2,ssh2);
                        % [~,latx,depth3,u3,dudt3,ududx3,vdudy3,wdudz3,zpgf3,corf3,frich3,fricv3] = zmb_isopycnal(lon,lat,depth,u3,v3,w3,temp3,salt3,ssh3);

                        if get_mean == 2

                            save MAM_quarter1_8000 lon latx depth1 temp1 u1 ssh1 taux1... 
                                dudt1 ududx1 vdudy1 wdudz1 zpgf1 corf1 frich1 fricv1

                            save MAM_quarter4_8000 lon latx depth4 temp4 u4 ssh4 taux4... 
                                dudt4 ududx4 vdudy4 wdudz4 zpgf4 corf4 frich4 fricv4
                            % 
                            % save MAM_quarter2 lon latx depth2 u2 ssh2 taux2... 
                            %     dudt2 ududx2 vdudy2 wdudz2 zpgf2 corf2 frich2 fricv2
                            % 
                            % save MAM_quarter3 lon latx depth3 u3 ssh3 taux3... 
                            %     dudt3 ududx3 vdudy3 wdudz3 zpgf3 corf3 frich3 fricv3
                        elseif get_mean == 3

                            save JJA_quarter1_8000 lon latx depth1 temp1 u1 ssh1 taux1... 
                                dudt1 ududx1 vdudy1 wdudz1 zpgf1 corf1 frich1 fricv1

                            save JJA_quarter4_8000 lon latx depth4 temp4 u4 ssh4 taux4... 
                                dudt4 ududx4 vdudy4 wdudz4 zpgf4 corf4 frich4 fricv4
                            %
                            % save JJA_quarter2 lon latx depth2 u2 ssh2 taux2... 
                            %     dudt2 ududx2 vdudy2 wdudz2 zpgf2 corf2 frich2 fricv2
                            % 
                            % save JJA_quarter3 lon latx depth3 u3 ssh3 taux3... 
                            %     dudt3 ududx3 vdudy3 wdudz3 zpgf3 corf3 frich3 fricv3
                            % 

                        end
                    end
end


