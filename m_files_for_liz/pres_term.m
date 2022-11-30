
function pres_term(salt_file, salt_var, temp_file, temp_var, ssh_file, ssh_var)

global dx dz rho0 g depth2 depth I_lon I_lat ntime nyearsub ncycle 

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
%     %
% salt = ncread(salt_file,salt_var, [I_lon(1) I_lat(1) 1 I_time(1)],...
%     [length(I_lon) length(I_lat) inf length(I_time)]);
% salt(salt<-10000) = NaN;
% 
% salt2 = interp1(depth,permute(salt,[3 2 1 4]), depth2);
% salt3 = permute([mean(salt2(:,1:2,:,:),2) mean(salt2(:,2:3,:,:),2) mean(salt2(:,3:4,:,:),2)],[3 2 1 4]);
% 
% salt=salt3; clear salt2 salt3
% 
% %
% temp = ncread(temp_file,temp_var,...
%     [I_lon(1) I_lat(1) 1 I_time(1)],[length(I_lon) length(I_lat) inf length(I_time)]);
% temp(temp<-10000) = NaN;
% 
% temp2 = interp1(depth, permute(temp,[3 2 1 4]), depth2);
% temp3 = permute([mean(temp2(:,1:2,:,:),2) mean(temp2(:,2:3,:,:),2) mean(temp2(:,3:4,:,:),2)],[3 2 1 4]);
% 
% temp=temp3; clear temp2 temp3
% 
% 
% % lon x lat x depth x time
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%
% % compute density 
% 
% salt2 = reshape(permute(salt, [3 1 2 4]),length(depth2),[]);
% clear salt
% 
% temp2 = reshape(permute(temp, [3 1 2 4]),length(depth2),[]);
% clear temp
% 
% rho = sw_dens(salt2, temp2, depth2);
% % 
% clear salt2 temp2
% 
% if j == 1
%    save rho rho -v7.3; 
%    
% else
%     rho_new = rho; load rho.mat; rho = cat(2,rho,rho_new);
%     save rho rho -v7.3; clear rho rho_new 
% end
% 
% clear rho
% 
% end
% % 
% %%
% %    for j = 1:ncycle
% % %   for j = 1:5
% % 
% % %for j = 6:11
% % disp(j)
% %     % set time bounds for given subset 
% %     % Important because need array sizes to stay within MATLAB memory limits
% %     
% %     if j ==ncycle
% %         I_time = (j-1)*nyearsub*12+1:ntime;
% %     else
% %         I_time = (j-1)*nyearsub*12+1:j*nyearsub*12;
% %     end
% %     
% %     %
% % salt = ncread(salt_file,salt_var, [I_lon(1) (I_lat(1)-1) 1 I_time(1)],...
% %     [length(I_lon) (length(I_lat)+2) inf length(I_time)]);
% % salt(salt<-10000) = NaN;
% % 
% % salt2 = interp1(depth,permute(salt,[3 2 1 4]), depth2);
% % salt3 = permute([mean(salt2(:,1:2,:,:),2) mean(salt2(:,2:3,:,:),2) mean(salt2(:,3:4,:,:),2) ...
% %     mean(salt2(:,4:5,:,:),2) mean(salt2(:,5:6,:,:),2)],[3 2 1 4] );
% % salt=salt3; clear salt2 salt3
% % 
% % %
% % temp = ncread(temp_file,temp_var,...
% %     [I_lon(1) (I_lat(1)-1) 1 I_time(1)],[length(I_lon) (length(I_lat)+2) inf length(I_time)]);
% % temp(temp<-10000) = NaN;
% % 
% % temp2 = interp1(depth, permute(temp,[3 2 1 4]), depth2);
% % temp3 = permute([mean(temp2(:,1:2,:,:),2) mean(temp2(:,2:3,:,:),2) mean(temp2(:,3:4,:,:),2) ...
% %     mean(temp2(:,4:5,:,:),2) mean(temp2(:,5:6,:,:),2)],[3 2 1 4]);
% % 
% % temp=temp3; clear temp2 temp3 I_time
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%
% % % compute potential density
% % 
% % salt2 = reshape(permute(salt, [3 1 2 4]),length(depth2),[]);
% % temp2 = reshape(permute(temp, [3 1 2 4]),length(depth2),[]);
% % 
% % clear salt temp
% % 
% % %rhoth = sw_pden(salt2,temp2, depth2, 0);
% % rhoth2 = sw_pden(salt2,temp2, depth2, 0);
% % 
% % clear salt2 temp2 
% % 
% % %     if j == 1
% % %        save rhoth rhoth;
% % %         else
% % %         rhoth_new = rhoth; load rhoth.mat; rhoth = cat(2,rhoth,rhoth_new);
% % %         save rhoth rhoth -v7; clear rhoth rhoth_new 
% % %         end
% %         
% %     if j == 6
% %        save rhoth2 rhoth2; 
% %        
% %     else
% %         rhoth2_new = rhoth2; load rhoth2.mat; rhoth2 = cat(2,rhoth2,rhoth2_new);
% %         save rhoth2 rhoth2 -v7; clear rhoth2 rhoth2_new 
% %     end
% % 
% %     % clear rhoth
% % clear rhoth2
% % 
% % end
% 
% 
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%
% for j = 1:ncycle
% disp(j)
%     % set time bounds for given subset 
%     % Important because need array sizes to stay within MATLAB memory limits
%     
%     k = nyearsub*12*3*length(I_lon);
%     
%     if j == ncycle
%         I_time_2 = (j-1)*k+1:ntime*3*length(I_lon);
%     else
%         I_time_2 = (j-1)*k+1:j*k;
%     end
%     
%     
%     if j == ncycle
%         I_time = (j-1)*nyearsub*12+1:ntime;
%     else
%         I_time = (j-1)*nyearsub*12+1:j*nyearsub*12;
%     end
% 
%     
% ssh = ncread(ssh_file, ssh_var, [I_lon(1) I_lat(1) I_time(1)],...
%     [length(I_lon) length(I_lat) length(I_time)]);
% ssh(ssh<-10000) = NaN;
% 
% ssh2 = [mean(ssh(:,1:2,:),2) mean(ssh(:,2:3,:),2) mean(ssh(:,3:4,:),2)];
% 
% ssh = ssh2; clear ssh2
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%
% % % compute pressure
% % 
% ssh2 = reshape(ssh,1,[]);
% clear ssh
% 
% load rho.mat
% 
% %%%%%%%% FIX RHO
% rho = rho(:,I_time_2);
% 
% p = permute(reshape(-1*g*(cumsum(rho)*dz + repmat((ssh2.*rho(1,:)),[198,1])),...
%     [length(depth2) length(I_lon) 3 length(I_time)]),[2 3 1 4]);
% 
% clear rho ssh2
% 
% if j == 1
%    save p p; 
%    
% else
%     p_new = p; load p.mat; p = cat(4,p,p_new);
%     save p p; clear p p_new 
% end
% 
% end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%compute zonal presure gradent force
% 
for j = 1:ncycle
 disp(j)
    
    if j == ncycle
        I_time = (j-1)*nyearsub*12+1:ntime;
    else
        I_time = (j-1)*nyearsub*12+1:j*nyearsub*12;
    end
 
        load p.mat

        p = squeeze(p(:,2,:,I_time));
    
% 
 zpgf = (-1/rho0)*((p(3:end,:,:)-p(1:end-2,:,:))/(2*dx));
% 
 clear p; 
 
    if j == 1
        save zpgf zpgf -v7.3; 
    %    
    else
        zpgf_new = zpgf; load zpgf.mat; zpgf = cat(3,zpgf,zpgf_new);
        save zpgf zpgf -v7.3;  clear zpgf zpgf_new 
    end
 
end
   
end
