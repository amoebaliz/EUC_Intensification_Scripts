
function fric_terms(u_file, u_var, temp_file, temp_var)

global dx dy dz depth2 depth I_lon I_lat ntime nyearsub ncycle Ah
% Constant Av options

% Av=ones(size(u))*1.66*10^(-3); % Bryden & Brady (1985)
% Av=ones(size(u))*1*10^(-4); % MITgcm manual

% VERTICALLY VARYING COEFFICIENT OF VERTICAL EDDY VISCOSITY (Av) SEE QIAO
% AND WEISBERG (1997, JPO, FIG. 14) THERE ARE TWO OPTIONS: TIE INFLECTION
% POINT TO THE THERMOCLINE (MAX DT/DZ) OR TO THE EUC (MAX U)


% Av_sfc=45*10^(-4); Av_tcline=3*10^(-4); Av_deep=15*10^(-4);

% for j = 1:ncycle
% disp(j)
%     % set time bounds for given subset 
%     % Important because need array sizes to stay within MATLAB memory limits
%     
%     if j ==ncycle
%         I_time = (j-1)*nyearsub*12+1:ntime;
%     else
%         I_time = (j-1)*nyearsub*12+1:j*nyearsub*12;
%     end
%     
% temp = ncread(temp_file,temp_var,...
%     [I_lon(1) I_lat(1) 1 I_time(1)],[length(I_lon) length(I_lat) inf length(I_time)]);
% temp(temp<-10000) = NaN;
% 
% temp2 = interp1(depth, permute(temp,[3 2 1 4]), depth2);
% temp3 = squeeze(permute(mean(temp2(:,2:3,:,:),2),[3 2 1 4]));
% 
% temp=temp3; clear temp2 temp3
% 
% Av=NaN(size(temp));
% 
% 
% % OPTION 1: TIE INFLECTION POINT TO THE THERMOCLINE
% 
% % deep is depth of thermocline + 100 m. note that QW97 says the values
% % below the EUC might be at least an order of magnitude smaller than this
% % based on measurements.
% 
% %
% for t = 1:length(I_time)
%     for x=1:size(Av,1)
%             % Identify the thermocline
%             z_tcline=min(depth2(diff(temp(x,:,t)) == min(diff(temp(x,:,t)))));
%             if( isfinite(z_tcline) && z_tcline ~= depth2(1))
%                 z_deep=z_tcline+100;
%                 a=[depth2(1) z_tcline z_deep];
%                 b=[Av_sfc Av_tcline Av_deep];
%                 c=depth2(1:find(depth2==z_deep));
%                 % d=spline(a,b,c);
%                 d=pchip(a,b,c);
%                 Av(x,1:find(depth2==z_deep),t)=d;
%                 Av(x,depth2>z_deep,t)=Av_deep;
%             end
%     end
% end; clear x y t z_tcline z_deep a b c d temp
% 
% 
% if j == 1
%     save Av Av
% else
%     Av_new = Av; load Av.mat
%     Av = cat(3,Av,Av_new); save Av Av
%     clear lat_inf Av Av_new 
% end
% end
%% Isopycnal conversions
% 

isopyc_layer_bnds = 1020:0.2:1028;


for j = 1:ncycle
disp(j)
    % set time bounds for given subset 
    % Important because need array sizes to stay within MATLAB memory limits
    
    if j ==ncycle
        I_time = (j-1)*nyearsub*12+1:ntime;
    else
        I_time = (j-1)*nyearsub*12+1:j*nyearsub*12;
    end

k = nyearsub*12*5*length(I_lon);
    
    if j == ncycle
        I_time_2 = (j-1)*k+1:ntime*5*length(I_lon);
    else
        I_time_2 = (j-1)*k+1:j*k;
    end
    
    
    %load Av.mat
    %Av = Av(:,:,I_time);
    
        if j <= 5
            load rhoth.mat
            I_time_rho = I_time;
            I_time_2_rho = I_time_2;
            
        else load rhoth2.mat
            rhoth = rhoth2;
            
            I_time_2_rho = I_time_2 - 5*k;
            I_time_rho = I_time - 5*nyearsub*12;
        end
        
    rhoth = rhoth(:,I_time_2_rho);
    rhoth = reshape(rhoth,length(depth2),length(I_lon), 5, length(I_time_rho));
    % rhoth = rhoth(:,:,3,:);
    rhoth = permute(rhoth(:,:,:,:),[2 3 1 4]);
    
    
u = ncread(u_file,u_var,[I_lon(1) (I_lat(1)-1) 1 I_time(1)],...
    [length(I_lon) (length(I_lat)+2) inf length(I_time)]);
u(u<-10000) = NaN;

u2 = interp1(depth,permute(u,[3 2 1 4]), depth2);
u3 = permute([mean(u2(:,1:2,:,:),2) mean(u2(:,2:3,:,:),2) mean(u2(:,3:4,:,:),2) ...
    mean(u2(:,4:5,:,:),2) mean(u2(:,5:6,:,:),2)],[3 2 1 4]);

u=u3; clear u2 u3

    % depth_iso = NaN(size(u,1),length(isopyc_layer_bnds)-1,size(u,4));
    % Av_iso=NaN(size(Av,1),length(isopyc_layer_bnds)-1,size(Av,3));
     u_iso=NaN(size(u,1),size(u,2),length(isopyc_layer_bnds)-1,size(u,4));
    
    %depth_array9{length(I_lon),length(isopyc_layer_bnds)-1,length(I_time)} = [];
    
    for x = 1:length(I_lon) % lon
        for y = 1:size(u,2) % lat
            for t = 1:length(I_time) % time
                for layer = 1: length(isopyc_layer_bnds)-1
                    zrange=depth2(squeeze(rhoth(x,y,:,t))>=isopyc_layer_bnds(layer)&squeeze(rhoth(x,y,:,t))<=isopyc_layer_bnds(layer+1));

                    if isfinite(zrange)
                         u_iso(x,y,layer,t)=nanmean(u(x,y,depth2>=zrange(1)&depth2<=zrange(length(zrange)),t),3);
                        
                         % if y == 3
                            % Av_iso(x,layer,t)=nanmean(Av(x,depth2>=zrange(1)&depth2<=zrange(length(zrange)),t),2);
                            % depth_iso(x,layer,t)=nanmean(depth2(depth2>=zrange(1)&depth2<=zrange(length(zrange))));
                            % depth_array9(x,layer,t) = {find(depth2>=zrange(1) & depth2<=zrange(end))};
                         % end
                        
                    end
                end
            end
        end
    end; clear x y xlayer rhoth zrange

u = u_iso; clear u_iso
% Av = Av_iso; clear Av_iso

dudx = squeeze((u(3:end,2,:,:)- u(1:end-2,2,:,:))/2*dx);
fricx = Ah*(dudx(3:end,:,:)-dudx(1:end-2,:,:)/2*dx);

clear dudx

dudy = (u(:,3:end,:,:)- u(:,1:end-2,:,:)/2*dy);
fricy = squeeze(Ah*(dudy(:,3,:,:)-dudy(:,1,:,:))/2*dy);

clear dudy

fricx2 = NaN(size(fricy)); fricx2(3:end-2,:,:)= fricx;
frich = fricx2 + fricy;

%% 

% save depth_array9 depth_array9

    


% 
%     if j == 1
%         save depth_array depth_array
%     else
%         depth_array_new = depth_array; load depth_array.mat; depth_array = cat(3,depth_array, depth_array_new);
%         save depth_array depth_array; clear depth_array depth_array_new
%     end
%     
%%
    if j == 1
        save frich frich
    else
        frich_new = frich; load frich.mat; frich = cat(3,frich, frich_new);
        save frich frich; clear frich frich_new
    end
    
%  dudz = squeeze(u(:,3,1:end-2,:)-u(:,3,3:end,:))./(depth_iso(:,3:end,:)-depth_iso(:,1:end-2,:));
%  Av_dudz = Av(:,2:end-1,:).*dudz; clear dudz
%  
%  fricv = (Av_dudz(:,1:end-2,:)-Av_dudz(:,3:end,:))./(depth_iso(:,4:end-1,:)-depth_iso(:,2:end-3,:));
%  clear Av_dudz
%  
% %%

%     if j == 1
%         save fricv fricv
%     else
%         fricv_new = fricv; load fricv.mat; fricv = cat(3,fricv, fricv_new);
%         save fricv fricv; clear fricv fricv_new
%     end
% %%%%%%%%%%%%%%%%%%%%%%%%% CHANGED TO VERT_FRIC.MAT 
%     
%     if j == 1
%         save depth_iso depth_iso
%     else
%         depth_iso_new = depth_iso; load depth_iso.mat; depth_iso = cat(3,depth_iso, depth_iso_new);
%         save depth_iso depth_iso; clear depth_iso depth_iso_new
%     end
%  

end
end
