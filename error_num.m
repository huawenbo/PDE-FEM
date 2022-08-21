%6点gauss数值积分计算误差函数
function [e0,e_1] = error_num(e,hx,hy,coord,connect,s)
S = hx*hy/2;
nodes = connect(e,:);%自由度编号
xe = coord(nodes,:); % 单元自由节点坐标
c1 = xe(2,1)*xe(3,2)-xe(3,1)*xe(2,2);
c2 = xe(3,1)*xe(1,2)-xe(1,1)*xe(3,2);
c3 = xe(1,1)*xe(2,2)-xe(2,1)*xe(1,2);
x1 = xe(2,2) - xe(3,2);x2 = xe(3,2) - xe(1,2);x3 = xe(1,2) - xe(2,2);
y1 = xe(3,1) - xe(2,1);y2 = xe(1,1) - xe(3,1);y3 = xe(2,1) - xe(1,1);
L = @(x,y) s(nodes(1))/(2*S)*(c1+x1*x+y1*y)+s(nodes(2))/(2*S)*(c2+x2*x+y2*y)+s(nodes(3))/(2*S)*(c3+x3*x+y3*y);
Lx=@(x,y) s(nodes(1))/(2*S)*x1+s(nodes(2))/(2*S)*x2+s(nodes(3))/(2*S)*x3;
Ly=@(x,y) s(nodes(1))/(2*S)*y1+s(nodes(2))/(2*S)*y2+s(nodes(3))/(2*S)*y3;
f=@(x,y) exp(x+y);
fx=@(x,y) exp(x+y);
fy=@(x,y) exp(x+y);
[x,y] = gauss(e,coord,connect);
f1=0;f2=0;
for i=1:6
    f1 = f1+S/6*(f(x(i),y(i))-L(x(i),y(i)))^2;
    f2 = f2+S/6*((fx(x(i),y(i))-Lx(x(i),y(i)))^2+(fy(x(i),y(i))-Ly(x(i),y(i)))^2);
end
e0 = f1;
e_1 = f2;
end