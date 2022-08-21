clc;%��������д���
clear;%��������ռ�
close all;%�ر�����ͼ��
%��������
tic
h=0.25;
k = 2/h;
nel = 3;%ÿ����Ԫ�Ľڵ���Ŀ,��ÿ����Ԫ���м����κ����������ã���Ԫ���ɶ�
Lx = 1;%���嵥Ԫ���ұ߽�
Ly = 1;%���嵥Ԫ���±߽�
N = k;%�ָ��һ������ĵ�Ԫ��Ŀ
paint(N);
numelx = N;%����ָ��x����Ԫ��Ŀ
numely = N;%����ָ��y����Ԫ��Ŀ
hx = 2*Lx/numelx;%x�����ϵĵ�Ԫ����
hy = 2*Ly/numely;%y�����ϵĵ�Ԫ����
numel = numelx*numely*2;%С��Ԫ����Ŀ,ÿ�����ηֳ�����������
numnodx = numelx + 1;%x����ڵ����
numnody = numely + 1;%y����ڵ����
numnod = numnodx*numnody;%�ܽڵ����
coordx = linspace(-Lx,Lx,numnodx)';%�ȷֽڵ������
coordy = linspace(-Ly,Ly,numnody)';%�ȷֽڵ������
[X,Y] = meshgrid(coordx,coordy);%�ų�����X��Y�ֱ��ʾ��Ӧλ�õĺ�������
X = X';Y = Y';coord = [X(:) Y(:)];%������չ��һ�С�
connect = connect_mat1(numnodx,numnody,nel);%���Ӿ��󣬱�ʾÿ����Ԫ��Χ�Ľڵ��ţ�Ҳ�����漰���κ�����ţ���ʱ��
ebcdof = unique([1 numnodx numnodx*numely+1:numnodx*numely+numnodx numnodx+1:numnodx:numnodx*(numely-1)+1 2*numnodx:numnodx:numely*numnodx]); % ǿ���Ա߽��ı��
bigk = sparse(numnod,numnod); %�նȾ���K
fext = sparse(numnod,1);%�غ�����f
%����ϵ������K���Ҷ���f��������λ�նȾ�����װ
for e = 1:numel %ͬһά���������Ȼ����Ԫ��ɨ��
    ke = elemstiff2d(e,nel,hx,hy,coord,connect);%���㵥Ԫ�նȾ���
    fe = elemforce2d(e,nel,hx,hy,coord,connect);%���㵥Ԫ�غ�����
    sctr = connect(e,:);
    bigk(sctr,sctr) = bigk(sctr,sctr) + ke;
    fext(sctr) = fext(sctr) + fe;
end
u_b=[];
for i=1:length(ebcdof)
    u_b=[u_b,exp(sum(coord(ebcdof(i),:)))];
end
ebcval = u_b; %����߽�ֵ��Ϊu_b
bound2=[2:2:2*(numnodx-1)];
for i=1:numnodx-1
    S = hx*hy/2;
    a=bound2(i);
    nodes=connect(a,:);
    xe = coord(nodes,:); % ��Ԫ���ɽڵ�����
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
%��ֵ��������
for i = 1:length(ebcdof)
    n = ebcdof(i);
    bigk(n,:) = 0;
    bigk(n,n) = 1;
    fext(n) = ebcval(i);
end
%�����ݶȷ���ⷽ��
u_coeff=bigk\fext;
% u_coeff = cg(bigk,fext,zeros(numnod,1),1e-50);%���ϵ������ʵ��Ҳ�Ǻ����ڶ�Ӧ���ϵ�ֵ
u_cal = u_coeff;
u_cal_re = reshape(u_coeff,numnodx,numnody);
%��ȷ��
L = Lx;
nsamp = 1001;
xsamp = linspace(-L,L,nsamp);%100�ȷ������м���100����
[X,Y] = meshgrid(xsamp,xsamp);
uexact = exactsolution2d(X(:),Y(:));
uexact_re = reshape(uexact,nsamp,nsamp);
%���ͼ
figure
h1=mesh(xsamp,xsamp,uexact_re);
title('Real Solutions');
xlabel('x');
ylabel('y');
% saveas(h1,['./plots/real solutions' '.jpg'])
%����Ԫ��ͼ
figure
h2 = mesh(coordx,coordy,u_cal_re);
title('FEM1 Solutions');
xlabel('x');
ylabel('y');
% saveas(h2,['./plots/' num2str(k) 'FEM1' '.jpg'])
%�������
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
