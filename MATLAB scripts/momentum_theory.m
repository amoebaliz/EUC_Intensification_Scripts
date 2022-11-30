close all
clear all
clc 

% time (e.g. month)
x = (1/12):(1/12):100;

% measurement (e.g. zonal velocity) 
y = sin(x)+x/3;

% derrivative (e.g. acceleration)
y2 = cos(x)+1/3;

% average acceleration NOTE THIS IS WHAT WE WERE DERIVING WITH THE MOMENTUM
% BALANCE!!!
trend_y = mean(y2);

% as a line: y = m*x
y_estimate = trend_y*x;

plot(x,y)
pause
hold on
plot(x,y_estimate,'k')
pause
% compare w/ trend_stat result
[slope2,plusminus,sig]=trend_stat(y',99);

% trend line: y= m*x
y_2estimate = slope2*(1:length(x));

plot(x,y_2estimate,'r')
