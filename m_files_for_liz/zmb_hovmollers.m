close all; clear all; clc

load ududx.mat; load ududx_trend.mat

nt = size(ududx,4); nmonth = 1:12:(nt-11); I_dep = find(depth2 == 51);

% initialize terms
ududx_hov_trend = NaN(12,length(slon(2:end-1))); ududx_hov_sig = ududx_hov_trend;

for m = 1:12
    I_month = (nmonth)+(m-1);
    for j = 1:length(slon(2:end-1))
        
        [ududx_hov_trend(m,j), ~,ududx_hov_sig(m,j)] = trend(squeeze(ududx(j,2,I_dep,I_month)),99);

    end
end

ududx_hov_trend(ududx_hov_sig==0)=NaN;
load redblue
colormap(redblue)

subplot(3,2,1)

[c,h] = contourf(slon(2:end-1),1:12,100*ududx_hov_trend,500);
set(h, 'EdgeColor','none')

hold on


cmax = 3*10^(-7);
cmin = -1*cmax;

caxis([cmin cmax])

colorbar
title('ududx')

clear all
%%
load vdudy.mat; load ududx_trend.mat

nt = size(vdudy,4); nmonth = 1:12:(nt-11); I_dep = find(depth2 == 51);

% initialize terms
vdudy_hov_trend = NaN(12,length(slon)); vdudy_hov_sig = vdudy_hov_trend;

for m = 1:12
    I_month = (nmonth)+(m-1);
    for j = 1:length(slon(2:end-1))
        
        [vdudy_hov_trend(m,j), ~,vdudy_hov_sig(m,j)] = trend(squeeze(vdudy(j,:,I_dep,I_month)),99);        
        
    end
end

vdudy_hov_trend(vdudy_hov_sig==0)=NaN;
load redblue
colormap(redblue)

subplot(3,2,3)

[c,h] = contourf(slon,1:12,100*vdudy_hov_trend,500);
set(h, 'EdgeColor','none')

hold on


cmax = 3*10^(-7);
cmin = -1*cmax;

caxis([cmin cmax])

colorbar

title('vdudy')

clear all
%%
load wdudz.mat; load ududx_trend.mat

nt = size(wdudz,4); nmonth = 1:12:(nt-11); I_dep = find(depth2 == 51);

wdudz_hov_trend = NaN(12,length(slon)); wdudz_hov_sig = wdudz_hov_trend;

for m = 1:12
    I_month = (nmonth)+(m-1);
    for j = 1:length(slon(2:end-1))
        
        [wdudz_hov_trend(m,j), ~,wdudz_hov_sig(m,j)] = trend(squeeze(wdudz(j,2,I_dep-1,I_month)),99);
                
    end
end

wdudz_hov_trend(wdudz_hov_sig==0)=NaN;
load redblue
colormap(redblue)

subplot(3,2,5)

[c,h] = contourf(slon,1:12,100*wdudz_hov_trend,500);
set(h, 'EdgeColor','none')

hold on


cmax = 3*10^(-7);
cmin = -1*cmax;

caxis([cmin cmax])

colorbar

title('wdudz')

clear all


clear all
%%
load zpgf.mat; load ududx_trend.mat

nt = size(zpgf,3); nmonth = 1:12:(nt-11); I_dep = find(depth2 == 51);

zpgf_hov_trend = NaN(12,length(slon(2:end-1))); zpgf_hov_sig = zpgf_hov_trend;

for m = 1:12
    I_month = (nmonth)+(m-1);
    for j = 1:length(slon(2:end-1))

        [zpgf_hov_trend(m,j), ~,zpgf_hov_sig(m,j)] = trend(squeeze(zpgf(j,I_dep,I_month)),99);
                
    end
end

zpgf_hov_trend(zpgf_hov_sig==0)=NaN;
load redblue
colormap(redblue)

subplot(3,2,2)

[c,h] = contourf(slon(2:end-1),1:12,100*zpgf_hov_trend,500);
set(h, 'EdgeColor','none')

hold on


cmax = 3*10^(-7);
cmin = -1*cmax;

caxis([cmin cmax])

colorbar

title('zpgf')
%%
clear all
%%
load frich.mat; load ududx_trend.mat

nt = size(frich,3); nmonth = 1:12:(nt-11); I_dep = 20;

frich_hov_trend = NaN(12,length(slon)); frich_hov_sig = frich_hov_trend;

for m = 1:12
    
    I_month = (nmonth)+(m-1);
    for j = 1:length(slon(2:end-1))

        [frich_hov_trend(m,j), ~,frich_hov_sig(m,j)] = trend(squeeze(frich(j,I_dep-2,I_month)),99);
        
    end
end

frich_hov_trend(frich_hov_sig==0)=NaN;
load redblue
colormap(redblue)

subplot(3,2,6)

[c,h] = contourf(slon,1:12,100*frich_hov_trend,500);
set(h, 'EdgeColor','none')

hold on


cmax = 3*10^(-6);
cmin = -1*cmax;

caxis([cmin cmax])

colorbar

title('frich')

%%
clear all

%%
clc
load fricv.mat; load ududx_trend.mat

nt = size(fricv,3); nmonth = 1:12:(nt-11); I_dep = 20;

fricv_hov_trend = NaN(12,length(slon)); fricv_hov_sig = fricv_hov_trend;

for m = 1:12
    I_month = (nmonth)+(m-1);
    for j = 1:length(slon)
        [fricv_hov_trend(m,j), ~,fricv_hov_sig(m,j)] = trend(squeeze(fricv(j,I_dep,I_month)),99);
        
    end
end

fricv_hov_trend(fricv_hov_sig==0)=NaN;
load redblue
colormap(redblue)

subplot(3,2,4)

[c,h] = contourf(slon,1:12,100*fricv_hov_trend,500);
set(h, 'EdgeColor','none')

hold on


cmax = 3*10^(-6);
cmin = -1*cmax;

caxis([cmin cmax])

colorbar

title('fricv')


%%
clear all

