%三角形元对应节点函数
function connect_mat = connect_mat1( numnodx,numnody,nel)
%输入横纵坐标的节点数目，和单元自由度
%输出连接矩阵，每个单元涉及的节点的编号
xn = 1:(numnodx*numnody);%拉成一条编号
A = reshape(xn,numnodx,numnody);%同形状编号
for j=1:numnody-1
    for i=1:numnodx-1
        a = A(i:i+1,j:j+1);
        a_vec = a(:);
        connect_mat(2*i+2*(j-1)*(numnody-1)-1:2*i+2*(j-1)*(numnody-1),1:nel) = [a_vec([1 4 3])';a_vec([4 1 2])'];
    end
end
end