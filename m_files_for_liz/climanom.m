function [Mc,Ma]=climanom(M)

% [CLIM,ANOM]=CLIMANOM(MONTHLY_DATA);
%
% Computes mean climatology and anomalies
%
% M:  Input monthly data matrix <time x lat x lon>
%
% Mc: Mean climatology <12 x lat x lon>
% Ma: Anomalies <time x lat x lon>
%
% For output climatology (Mc) to be ordered Jan -> Dec, first month of
% input data matrix (M) must be Jan!

% compute monthly climatology
nt=size(M,1);
nyears=floor(nt/12);
nxmo=nt-nyears*12;
%nxmo = rem(nt,12);
Mc=zeros(12,size(M,2),size(M,3));
for icalmo=1:12
    Mc(icalmo,:,:)=nanmean(M([icalmo:12:nt],:,:),1);
end

% remove monthly climatology for monthly anomalies
 % Ma = M -[repmat(Mc,[nyears 1 1]); Mc(1:rem(nt,12),:,:)];

Mc_rep=zeros(nt,size(M,2),size(M,3));
for iyear=1:nyears
    Mc_rep((iyear-1)*12+1:(iyear-1)*12+12,:,:)=Mc(:,:,:);
end

Ma=zeros(nt,size(M,2),size(M,3));
for t=1:nt
    Ma(t,:,:)=M(t,:,:)-Mc_rep(t,:,:);
end