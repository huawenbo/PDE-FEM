%共轭梯度法
function [x,k]=cg(A,b,x0,delta)
h=size(A,1);
x=zeros(h,1);
r=zeros(h,1);
p=zeros(h,1);
a=[];
c=[];
x(:,1)=x0;
r(:,1)=b-A*x0;
k=1;
while sum(r(:,k)==zeros(h,1))~=h && r(:,k)'*r(:,k)>delta %终止条件。
    k=k+1;
    if k==2
        p(:,1)=r(:,1);
    else
        c(k-2)=r(:,k-1)'*r(:,k-1)/(r(:,k-2)'*r(:,k-2));
        p(:,k-1)=r(:,k-1)+c(k-2)*p(:,k-2);
    end
    a(k-1)=r(:,k-1)'*r(:,k-1)/(p(:,k-1)'*A*p(:,k-1));
    x(:,k)=x(:,k-1)+a(k-1)*p(:,k-1);
    r(:,k)=r(:,k-1)-a(k-1)*A*p(:,k-1);
end
x=x(:,k);
k=k-1;
end