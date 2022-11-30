function [trend,plusminus,sig,trend_pom]=trend_stat(M,siglev)

% [trend,plusminus,sig]=trend(M,siglev);
% 
% INPUT:
% M is the input time series <nt x 1>
% siglev is the signficance level (e.g., 95)
% 
% OUTPUT:
% trend is the linear trend
% plusminus is the value for trend +/- plusminus
% sig is 1 for significant, 0 for not significant
% trend_pom is trend expressed as percent of mean (pom)
%
% NOTES:
%
% Simply uses MATLAB's linear regression function to calculate a trend.
%
% Does not explicitly take into account degrees of freedom, so
% significance may not be valid here if time series is smoothed/filtered.
%
% Should still work if there are NaNs in the input time series.
%
% The units of trend are the units of the values in the input time series
% per unit increment of time between values in the input time series.
%
% Time spacing must be linear.

alpha=1-siglev/100;
X=ones(size(M,1),2);
X(:,2)=[1:size(M,1)]';

yM=M;
[b,bint]=regress(yM,X,alpha);
trend=b(2);

sig=0;
plusminus=(bint(2,2)-bint(2,1))/2;
if ( abs(trend) > abs(plusminus) )
    sig=1;
end

trend_pom=100*trend/mean(M);