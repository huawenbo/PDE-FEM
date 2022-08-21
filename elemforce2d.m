%载荷向量数值计算函数，采用6点gauss数值积分
function [fe] = elemforce2d(e,nel,hx,hy,coord,connect)
S = hx*hy/2;
fe = zeros(nel,1); %初始化载荷向量
nodes = connect(e,:);%自由度编号
xe = coord(nodes,:); % 单元自由节点坐标
c1 = xe(2,1)*xe(3,2)-xe(3,1)*xe(2,2);
c2 = xe(3,1)*xe(1,2)-xe(1,1)*xe(3,2);
c3 = xe(1,1)*xe(2,2)-xe(2,1)*xe(1,2);
x1 = xe(2,2) - xe(3,2);x2 = xe(3,2) - xe(1,2);x3 = xe(1,2) - xe(2,2);
y1 = xe(3,1) - xe(2,1);y2 = xe(1,1) - xe(3,1);y3 = xe(2,1) - xe(1,1);
L1 = @(x,y) 1/(2*S)*(c1+x1*x+y1*y);
L2 = @(x,y) 1/(2*S)*(c2+x2*x+y2*y);
L3 = @(x,y) 1/(2*S)*(c3+x3*x+y3*y);
[x,y] = gauss(e,coord,connect);
f1=0;f2=0;f3=0;
for i=1:6
    f1 = f1+S/6*fun(x(i),y(i))*L1(x(i),y(i));
    f2 = f2+S/6*fun(x(i),y(i))*L2(x(i),y(i));
    f3 = f3+S/6*fun(x(i),y(i))*L3(x(i),y(i));
end
fe = [f1;f2;f3];
end