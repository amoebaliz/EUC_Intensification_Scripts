function du_cor_terms(u_file, u_var, v_file, v_var, w_file, w_var)

global dx dy dz depth2 depth I_lon I_lat ntime nyearsub ncycle

for j = 1:ncycle

    % set time bounds for given subset 
    % Important because need array sizes to stay within MATLAB memory limits
    
    if j ==ncycle
        I_time = (j-1)*nyearsub*12+1:ntime;
    else
        I_time = (j-1)*nyearsub*12+1:j*nyearsub*12;
    end

% Get data from SODA files

% 
% u = ncread(u_file,u_var,[I_lon(1) I_lat(1) 1 I_time(1)],...
%     [length(I_lon) length(I_lat) inf length(I_time)]);
% u(u<-10000) = NaN;
% 
% u2 = interp1(depth,permute(u,[3 2 1 4]), depth2);
% u3 = permute([mean(u2(:,1:2,:,:),2) mean(u2(:,2:3,:,:),2) mean(u2(:,3:4,:,:),2)],[3 2 1 4]);
% 
% u=u3; clear u2 u3

%
v = ncread(v_file,v_var, [I_lon(1) I_lat(1) 1 I_time(1)],...
    [length(I_lon) length(I_lat) inf length(I_time)]);
v(v<-10000) = NaN;

v2 = interp1(depth,permute(v,[3 2 1 4]), depth2);
v3 = permute([mean(v2(:,1:2,:,:),2) mean(v2(:,2:3,:,:),2) mean(v2(:,3:4,:,:),2)],[3 2 1 4]);

v=v3; clear v2 v3

% %
% w = ncread(w_file,w_var,[I_lon(1) I_lat(1) 1 I_time(1)],...
%     [length(I_lon) length(I_lat) inf length(I_time)]);
% w(w<-10000) = NaN;
% 
% w2 = interp1(depth,permute(w,[3 2 1 4]), depth2);
% w3 = permute([mean(w2(:,1:2,:,:),2) mean(w2(:,2:3,:,:),2) mean(w2(:,3:4,:,:),2)],[3 2 1 4]);
% w=w3; clear w2 w3

% all are in lon lat (depth) time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
% 
% % compute nonlinear advective terms
% 
% ududx = u(2:end-1,:,:,:).*(u(3:end,:,:,:)-u(1:end-2,:,:,:))/(2*dx);
% 
% if size(v,2) == 3
%     vdudy = v(:,2,:,:).*(u(:,3,:,:)-u(:,1,:,:))/(2*dy);
% else
%     vdudy = v(:,2:end-1,:,:).*(v(:,3:end,:,:)-v(:,1:end-2,:,:))/(2*dy);
% end
% 
% wdudz = w(:,:,2:end-1,:).*(u(:,:,3:end,:)-u(:,:,1:end-2,:))/(2*dz);
% 
% if j == 1
%    save ududx ududx; save vdudy vdudy; save wdudz wdudz
%    clear ududx vdudy wdudz
%    
% else
%     ududx_new = ududx; load ududx.mat; ududx = cat(4,ududx,ududx_new);
%     save ududx ududx; clear ududx ududx_new 
%     
%     vdudy_new = vdudy; load vdudy.mat; vdudy = cat(4,vdudy,vdudy_new);
%     save vdudy vdudy; clear vdudy vdudy_new 
%     
%     wdudz_new = wdudz; load wdudz.mat; wdudz = cat(4,wdudz,wdudz_new);
%     save wdudz wdudz; clear wdudz wdudz_new 
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
% compute coriolis force
omega=7.292*10^(-5);

lat_inf = permute(repmat(sind(-.5:.5:.5)',[1 size(v,1) size(v,3) size(v,4)]),[2 1 3 4]);
corf = 2*omega*lat_inf.*v; 

if j == 1
    save corf corf
else
    corf_new = corf; load corf.mat
    corf = cat(4,corf,corf_new); save corf corf
    clear lat_inf corf corf_new 
 
end

end 