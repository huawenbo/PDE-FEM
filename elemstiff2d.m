%单位刚度矩阵数值计算函数
function [ke] = elemstiff2d(e,nel,hx,hy,coord,connect)
S = hx*hy/2;
ke = zeros(nel,nel); %单刚初始化
nodes = connect(e,:);%相关形函数编号
xe = coord(nodes,:); %相关节点的坐标
x1 = xe(2,2) - xe(3,2);x2 = xe(3,2) - xe(1,2);x3 = xe(1,2) - xe(2,2);
y1 = xe(3,1) - xe(2,1);y2 = xe(1,1) - xe(3,1);y3 = xe(2,1) - xe(1,1);
auv11 = (x1*x1 + y1*y1)/(4*S)+2*S*factorial(2)/factorial(4);
auv12 = (x1*x2 + y1*y2)/(4*S)+2*S*factorial(1)/factorial(4);
auv13 = (x1*x3 + y1*y3)/(4*S)+2*S*factorial(1)/factorial(4);
auv22 = (x2*x2 + y2*y2)/(4*S)+2*S*factorial(2)/factorial(4);
auv23 = (x2*x3 + y2*y3)/(4*S)+2*S*factorial(1)/factorial(4);
auv33 = (x3*x3 + y3*y3)/(4*S)+2*S*factorial(2)/factorial(4);
ke = ke + [auv11 auv12 auv13;
           auv12 auv22 auv23;
           auv13 auv23 auv33];
end