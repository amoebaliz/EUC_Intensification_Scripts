




load vert_fric.mat


load depth_array.mat
depth_array = depth_array(:,3:end-2,:);

depth2 = 5:2:400;
%%
% initialize new vairable
vert_fric_dep = NaN(size(fricv,1),length(depth2),size(fricv,3));
%%


    for j = 1:size(fricv,1)             % longitude only 223
        disp(j)
        for k = 1:size(fricv,2)         % depth - only 36 vs 40
            for t = 1:size(fricv,3)     % time - 1716 ... for now

                if isfinite(depth_array{j,k,t})

                    i_deps = cell2mat(depth_array(j,k,t));
                    % assign the average friction value to those depths 
                    vert_fric_dep(j,i_deps,t) = fricv(j,k,t);

                end

            end
        end
       
    end
    
    
    save vert_fric_dep vert_fric_dep
    
    
  
    