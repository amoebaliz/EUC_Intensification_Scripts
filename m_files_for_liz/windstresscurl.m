function crl=windstresscurl(lon,lat,taux,tauy)

% compute wind stress curl
% 
% NOTES:
% 1. taux and tauy must be <time x lat x lon>
% 2. if units of taux and tauy are N m^-2 (=Pa), units of crl are 10^-8 N m^-3
% 3. be careful that taux and tauy are positive eastward and northward, respectively. sometimes stress is expressed in the opposite sense as the wind direction.
%
% compare in the tropical Pacific to:
% http://faculty.washington.edu/kessler/pacs/inverse/sverers-9197.gif
% http://journals.ametsoc.org/doi/full/10.1175/1520-0469%282001%29058%3C0109%3ATGDOTT%3E2.0.CO%3B2

crl=zeros(size(taux));
for t=1:size(taux,1);
    [crl(t,:,:),~]=curl(lon,lat,squeeze(taux(t,:,:)),squeeze(tauy(t,:,:)));
end; clear t
crl=crl*10^3;