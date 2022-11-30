function F = filtrage2(X,type,N,per)

% F = FILTRAGE2(X,TYPE,N,PER)
%
% This function filters the columns of an input matrix with a recursive
% Butterworth filter.
%
% Input arguments...
%
% X:    A matrix of real numbers <time x station>.
%       The columns of X are independently filtered
% type: A character string to designate the filter type
%       low, high, bandpass, or stop
%       If bandpass or stop, per should be a two-member vector.
% N:    Odd scalar giving the order of the filter
% per:  A scalar or vector giving the periods to define the filter
%       For annual data and per=50, the cutoff frequency
%       is 1/50 cycles/year. For monthly data and per=120, the cutoff
%       frequency is 1/120 cycles/month or 1/10 cycles per year.
% 
% Outputs...
%
% F: Filtered matrix <time x station>
%
% created in May 1996 and modified in May 2004
% Vincent Moron
% moron@cerege.fr
%
% Modified by Kristopher B. Karnauskas
% Woods Hole Oceanographic Institution
% kk@whoi.edu - February 2011
%
% mostly just corrected notation/cleaned up code and doc
% preallocated bf

[n,nc]=size(X);
ns=round(n/2);
sim1=flipud(X(2:ns+1,:));
sim2=flipud(X(n-ns:n-1,:));
X=[sim1;X;sim2];

if strcmp(type,'high') | strcmp(type,'low');
    per=per(1);
    freq=2./per;
elseif strcmp(type,'stop') | strcmp(type,'bandpass');
    freq=2./per;
end

[b,a]=butter(N,freq,type);

bf=zeros(size(X));
for i=1:nc,
    bf(:,i)=filtfilt(b,a,X(:,i));
end;

F=bf(ns+1:n+ns,:);

