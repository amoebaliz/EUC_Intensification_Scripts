% year=linspace(1982,2011-(1/10613),10613);
% year_nino3=linspace(1982,2011-(1/12),348);

figure(1)
subplot(3,2,1)
    plot(year,dsst)
    xlim([1982 2011])
    ylim([-1.6 1.6])
    hold on
    trend1=dsst-detrend(dsst);
    plot(year,trend1,'c')
    title((trend1(2)-trend1(1))*365*100);
subplot(3,2,2)
    dsst2=runmean(dsst,31);
    plot(year,dsst2)
    xlim([1982 2011])
    ylim([-0.8 0.8])
    hold on
    year2=year;year2(isnan(dsst2))=[];
    dsst2(isnan(dsst2))=[];
    trend2=dsst2-detrend(dsst2);
    plot(year2,trend2,'c')
    title((trend2(2)-trend2(1))*365*100);
subplot(3,2,3)
    dsst3=runmean(dsst,365);
    plot(year,dsst3)
    xlim([1982 2011])
    ylim([-0.4 0.4])
    hold on
    year3=year;year3(isnan(dsst3))=[];
    dsst3(isnan(dsst3))=[];
    trend3=dsst3-detrend(dsst3);
    plot(year3,trend3,'c')
    title((trend3(2)-trend3(1))*365*100);
subplot(3,2,4)
    dsst4=runmean(dsst,365*7);
    plot(year,dsst4)
    xlim([1982 2011])
    ylim([-0.2 0.2])
    hold on
    year4=year;year4(isnan(dsst4))=[];
    dsst4(isnan(dsst4))=[];
    trend4=dsst4-detrend(dsst4);
    plot(year4,trend4,'c')
    title((trend4(2)-trend4(1))*365*100);
subplot(3,2,5)
    nino33=runmean(nino3,13);
    plot(year_nino3,nino33,'r')
    xlim([1982 2011])
    ylim([-3 3])
    hold on
    year_nino33=year_nino3;year_nino33(isnan(nino33))=[];
    nino33(isnan(nino33))=[];
    trend_nino33=nino33-detrend(nino33);
    plot(year_nino33,trend_nino33,'y')
    title((trend_nino33(2)-trend_nino33(1))*12*100);
    
    
    