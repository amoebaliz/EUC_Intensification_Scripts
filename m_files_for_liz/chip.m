close all; clear all
clc

tic
% identify the file 
fid = fopen('SST-STR_export_data_20130415.txt');
%fid = fopen('test.txt');


% get the data out of the file
% NOTE: you MUST specify what is in each column, whether text or numbers
%                  numbers = %n
%                  text    = %s

C = textscan(fid, '%s%s%n%n%n%s%s%s%s%n%s%s%s%s%s%s%s%s%s%s%s%s%n%n%n%n%n%n%n%n%n',700,'Delimiter',',');

%C = textscan(fid, '%n%n%s');

% close up the main file
fclose(fid);
toc

%% Assign Each Data Column to a Structure

NOAA_data.number = C{1};
NOAA_data.data = C{2};
NOAA_data.date = C{3};


%%

Index_values = 1:length(C{1});

%%%%%%%%
%%%%%%%%                                    |             |   SPECIFY DATA
%%%%%%%%  INPUT YOUR SEARCH TEXT HERE      _|_           _|_  COLUMN HERE
%%%%%%%%                                   \ /           \ /
%%%%%%%%                                    V             V

[i]=ind2sub(size(NOAA_data.date), strcmp('Thurs', NOAA_data.date));
[j]=ind2sub(size(NOAA_data.date), strcmp('Mon', NOAA_data.date));
[k]=ind2sub(size(NOAA_data.date), strcmp('Tues', NOAA_data.date));
[m]=ind2sub(size(NOAA_data.date), strcmp('Wed', NOAA_data.date));


% Get the numberical index values for the data
I_data1 = Index_values(i == 1 | j == 1);
I_data2 = Index_values(k == 1);
I_data3 = Index_values(m == 1);

% Use the indeces picked out by the text search to get the specific data
% you want to plot

XData1 = NOAA_data.number(I_data1);
YData1 = NOAA_data.data(I_data1); 

XData2 = NOAA_data.number(I_data2);
YData2 = NOAA_data.data(I_data2);

XData3 = NOAA_data.number(I_data3);
YData3 = NOAA_data.data(I_data3);

%%%%%%%%%             Plotting up the actual data              %%%%%%%%%%%%
%
% For info on adjusting plot styles: 
%           http://www.mathworks.com/help/matlab/ref/linespec.html

plot(XData1, YData1,'d','MarkerFaceColor','y','MarkerEdgeColor','k',...
    'MarkerSize',10)
hold on
plot(XData2, YData2,'bv','MarkerFaceColor','g','MarkerSize',10)
plot(XData3, YData3,'ms','MarkerSize',10)


% Messing with the plot labels
xlabel('X axis','FontSize', 14)
ylabel('Y axis','FontSize', 14)