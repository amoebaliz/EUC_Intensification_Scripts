close all; clear all; clc

load redblue

x = 0:.05:2*pi;
y = sin(x);

yrange = linspace(min(y), max(y),length(redblue));

for j = 1:length(x)
    
    rb_I1 = find(yrange<=y(j),1,'last');
    rb_I2 = find(yrange>=y(j),1,'first');
    
    interval = yrange(rb_I2)-yrange(rb_I1);
    
    if interval ~= 0
        frac1 = (y(j) -  yrange(rb_I1))/interval;
        frac2 = (-1*(y(j) - yrange(rb_I2)))/interval;

        colval = redblue(rb_I1,:)*(1-frac1) + redblue(rb_I2,:)*(1-frac2);
    
    else 
        colval = redblue(rb_I1,:);
    end
    
    colval_2 = colval/2;
        
    plot(x(j),y(j),'o', 'MarkerEdgeColor', colval_2,...
        'MarkerFaceColor', colval) 
    
    xlim([0 2*pi])
    ylim([-1 1])
    
    hold on
    pause(0.01)

end

colormap(redblue)
colorbar
caxis([-1, 1])

