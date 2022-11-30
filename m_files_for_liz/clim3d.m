function Mc=clim3d(M)

% [CLIM]=CLIM3D(MONTHLY_DATA);
%
% Computes mean climatology of 3D data (actually 4D including time)
%
% M:  Input monthly data matrix <time x lev x lat x lon>
%
% Mc: Mean climatology <12 x lev x lat x lon>
%
% For output climatology (Mc) to be ordered Jan -> Dec, first month of
% input data matrix (M) must be Jan!

% compute monthly climatology
nt=size(M,1);
nyears=floor(nt/12);
nxmo=nt-nyears*12;
Mc=zeros(12,size(M,2),size(M,3),size(M,4));
for icalmo=1:12
    Mc(icalmo,:,:,:)=nanmean(M([icalmo:12:nt],:,:,:),1);
end