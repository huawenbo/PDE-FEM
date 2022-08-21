%三角剖分图解函数
function paint(N)
syms x
x=-1:2:1;
y=-1:2:1;
figure
for h=1:N
    y1=x+2*(h-1)/N;
    y2=x-2*(h-1)/N;
    y3=[2*(h-N/2)/N,2*(h-N/2)/N];
    x4=[2*(h-N/2)/N,2*(h-N/2)/N];
    plot(x,y1,'black')
    hold on
    plot(x,y2,'black')
    hold on
    plot(x,y3,'black')
    hold on
    plot(x4,y,'black')
    hold on
end
title('Triangulation diagram')
xlabel('x')
ylabel('y')
axis([-1,1,-1,1]);
% saveas(gcf,['./plots/' num2str(N) 'Triangulation_diagram' '.jpg'])
end
