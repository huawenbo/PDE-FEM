clc;%清空命令行窗口
clear;%清除工作空间
close all;%关闭所有图像
%参数设置
tic
h=0.25;
k = 2/h;
nel = 3;%每个单元的节点数目,即每个单元上有几个形函数参与作用，单元自由度
Lx = 1;%定义单元左右边界
Ly = 1;%定义单元上下边界
N = k;%分割的一个方向的单元数目
paint(N);
numelx = N;%定义分割的x方向单元数目
numely = N;%定义分割的y方向单元数目
hx = 2*Lx/numelx;%x方向上的单元长度
hy = 2*Ly/numely;%y方向上的单元长度
numel = numelx*numely*2;%小单元的数目,每个矩形分成两个三角形
numnodx = numelx + 1;%x方向节点个数
numnody = numely + 1;%y方向节点个数
numnod = numnodx*numnody;%总节点个数
coordx = linspace(-Lx,Lx,numnodx)';%等分节点的坐标
coordy = linspace(-Ly,Ly,numnody)';%等分节点的坐标
[X,Y] = meshgrid(coordx,coordy);%张成网格，X和Y分别表示对应位置的横纵坐标
X = X';Y = Y';coord = [X(:) Y(:)];%将坐标展成一列。
connect = connect_mat1(numnodx,numnody,nel);%连接矩阵，表示每个单元周围的节点编号，也就是涉及的形函数编号，逆时针
ebcdof = unique([1 numnodx numnodx*numely+1:numnodx*numely+numnodx numnodx+1:numnodx:numnodx*(numely-1)+1 2*numnodx:numnodx:numely*numnodx]); % 强制性边界点的编号
bigk = sparse(numnod,numnod); %刚度矩阵K
fext = sparse(numnod,1);%载荷向量f
%计算系数矩阵K和右端项f，并将单位刚度矩阵组装
for e = 1:numel %同一维的情况，依然按单元来扫描
    ke = elemstiff2d(e,nel,hx,hy,coord,connect);%计算单元刚度矩阵
    fe = elemforce2d(e,nel,hx,hy,coord,connect);%计算单元载荷向量
    sctr = connect(e,:);
    bigk(sctr,sctr) = bigk(sctr,sctr) + ke;
    fext(sctr) = fext(sctr) + fe;
end
u_b=[];
for i=1:length(ebcdof)
    u_b=[u_b,exp(sum(coord(ebcdof(i),:)))];
end
ebcval = u_b; %假设边界值都为u_b
bound2=[2:2:2*(numnodx-1)];
for i=1:numnodx-1
    S = hx*hy/2;
    a=bound2(i);
    nodes=connect(a,:);
    xe = coord(nodes,:); % 单元自由节点坐标
    c1 = xe(2,1)*xe(3,2)-xe(3,1)*xe(2,2);
    c2 = xe(3,1)*xe(1,2)-xe(1,1)*xe(3,2);
    c3 = xe(1,1)*xe(2,2)-xe(2,1)*xe(1,2);
    x1 = xe(2,2) - xe(3,2);x2 = xe(3,2) - xe(1,2);x3 = xe(1,2) - xe(2,2);
    y1 = xe(3,1) - xe(2,1);y2 = xe(1,1) - xe(3,1);y3 = xe(2,1) - xe(1,1);
    L1 = @(x,y) 1/(2*S)*(c1+x1*x+y1*y);
    L2 = @(x,y) 1/(2*S)*(c2+x2*x+y2*y);
    L3 = @(x,y) 1/(2*S)*(c3+x3*x+y3*y);
    t=[0 -0.538469310106 0.538469310106 -0.906179845939 0.906179845939];
    A=[0.56888888889 0.478628670499 0.478628670499 0.236926885056 0.236926885056];
    F=@(x) -exp(x-1);
    x=(xe(3,1)-xe(2,1))/2*t+(xe(3,1)+xe(2,1))/2;
    f1=(L2(x,-1).*F(x))*A'*(xe(3,1)-xe(2,1))/2;
    f2=(L3(x,-1).*F(x))*A'*(xe(3,1)-xe(2,1))/2;
    fext(nodes(2)) = fext(nodes(2)) + f1;
    fext(nodes(3)) = fext(nodes(3)) + f2;
end
%边值条件处理
for i = 1:length(ebcdof)
    n = ebcdof(i);
    bigk(n,:) = 0;
    bigk(n,n) = 1;
    fext(n) = ebcval(i);
end
%共轭梯度法求解方程
u_coeff=bigk\fext;
% u_coeff = cg(bigk,fext,zeros(numnod,1),1e-50);%求出系数，事实上也是函数在对应点上的值
u_cal = u_coeff;
u_cal_re = reshape(u_coeff,numnodx,numnody);
%求精确解
L = Lx;
nsamp = 1001;
xsamp = linspace(-L,L,nsamp);%100等分区间中间有100个数
[X,Y] = meshgrid(xsamp,xsamp);
uexact = exactsolution2d(X(:),Y(:));
uexact_re = reshape(uexact,nsamp,nsamp);
%真解图
figure
h1=mesh(xsamp,xsamp,uexact_re);
title('Real Solutions');
xlabel('x');
ylabel('y');
% saveas(h1,['./plots/real_solutions' '.jpg'])
%有限元解图
figure
h2 = mesh(coordx,coordy,u_cal_re);
title('FEM1 Solutions');
xlabel('x');
ylabel('y');
% saveas(h2,['./plots/' num2str(k) 'FEM1' '.jpg'])
%计算误差
u_ex = exactsolution2d(coord(:,1),coord(:,2));
u_ex_re = reshape(u_ex,numnodx,numnody);
e=0;
e2=0;
for i=1:numel
    [a,b]=error_num(i,hx,hy,coord,connect,u_cal);
    e=e+a;
    e2=e2+b;
end
error_L2 = sqrt(e) 
error_diffL2 = sqrt(e2)
error_H0 = sqrt(e+e2)
toc
