%6点gauss数值积分节点函数
function [x,y]=gauss(e,coord,connect)
nodes = connect(e,:);%自由度编号
xe = coord(nodes,:); %单元自由节点坐标
L1=0.659027622374092;
L2=0.231933368553031;
L3=0.109039009072877;
x1=xe(1,1)*L1+xe(2,1)*L2+xe(3,1)*L3;
x2=xe(1,1)*L1+xe(2,1)*L3+xe(3,1)*L2;
x3=xe(1,1)*L2+xe(2,1)*L1+xe(3,1)*L3;
x4=xe(1,1)*L2+xe(2,1)*L3+xe(3,1)*L1;
x5=xe(1,1)*L3+xe(2,1)*L1+xe(3,1)*L2;
x6=xe(1,1)*L3+xe(2,1)*L2+xe(3,1)*L1;
y1=xe(1,2)*L1+xe(2,2)*L2+xe(3,2)*L3;
y2=xe(1,2)*L1+xe(2,2)*L3+xe(3,2)*L2;
y3=xe(1,2)*L2+xe(2,2)*L1+xe(3,2)*L3;
y4=xe(1,2)*L2+xe(2,2)*L3+xe(3,2)*L1;
y5=xe(1,2)*L3+xe(2,2)*L1+xe(3,2)*L2;
y6=xe(1,2)*L3+xe(2,2)*L2+xe(3,2)*L1;
x=[x1,x2,x3,x4,x5,x6];
y=[y1,y2,y3,y4,y5,y6];
end