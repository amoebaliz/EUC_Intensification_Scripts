function [xcoast,ycoast] = getcoast

load coast
long(long>=-180 & long<0)=long(long>=-180 & long<0)+360;
lat(long>359.5)=NaN;lat(long<0.5)=NaN;
long(long>359.5)=NaN;long(long<0.5)=NaN;

xcoast=long;
ycoast=lat;