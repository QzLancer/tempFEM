x = 1:0.01:2.5;
plot(x,sin(x)-x/2,'Color','black');
hold on;
line([1,2.5],[0,0],'color','red');
hold on;
x0 = 1.4;
y0 = sin(x0)-x0/2;
scatter(x0,y0,'filled','black');
hold on;

for i = 1:20
    x1 = x0 - (sin(x0)-x0/2)/(cos(x0)-1/2);
    y1 = sin(x1)-x1/2;
    scatter(x1,y1,'filled','black');
    hold on;
    line([x0,x1],[y0,0],'Color','black');
    hold on;
    line([x1,x1],[0,y1],'LineStyle','--','Color','black');
    x0 = x1;
    y0 = y1;
end