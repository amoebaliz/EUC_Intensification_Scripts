function Mc=clim2d(M)

% [CLIM]=CLIM2D(MONTHLY_DATA);
%
% Computes mean climatology of 2D data (actually 3D including time)
%
% M:  Input monthly data matrix <time x lat x lon>
%
% Mc: Mean climatology <12 x lat x lon>
%
% For output climatology (Mc) to be ordered Jan -> Dec, first month of
% input data matrix (M) must be Jan!

% compute monthly climatology
nt=size(M,1);
nyears=floor(nt/12);
nxmo=nt-nyears*12;
Mc=zeros(12,size(M,2),size(M,3));
for icalmo=1:12
    Mc(icalmo,:,:)=nanmean(M([icalmo:12:nt],:,:),1);
end
