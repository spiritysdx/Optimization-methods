A=[1 -2 1;4 -1 -2;-2 0 1];
B=[11 -3 1];
B=B';
C=[3 -1 -1];
C=C';
target=1;
sign=[-1 -1 0];
sign=sign';
[x,y,result]=TwoStageSimplexAlgorithm(A,B,C,target,sign);
function [x,y,ResultFlag]=TwoStageSimplexAlgorithm(A,B,C,target,sign)
% 2020-4-2 臻orz
%inputs:
%   A:系数矩阵 m*n
%   B:右端向量 m*1
%   C:价格系数向量 n*1
%alternative inputs:
%   target:优化目标 0 ~ min; 1 ~ max;
%   sign:约束条件符号向量 m*1 其中
%       -1 ~ '<=';0 ~ '=';1 ~ '>'; 默认为0,即为等式
%outputs:
%   x:最优解 n*1
%   y:最优值 num
%   ResultFlag:是否找到最优解
[m,n] = size(A);
P = [];
x = zeros(n,1);
y = 0;
ResultFlag = 0;
j = 0;

%standardization
if target
    C = -C;%目标函数的转化
end
A(B<0,:) = -A(B<0,:);
sign(B<0,:) = -sign(B<0,:);
B = abs(B);%约束条件的转化
for i = sign'
    j = j+1;
    switch i
        case -1%引入松弛变量
            a = zeros(m,1);a(j) = 1;
            A = [A a];
            C = [C;0];
        case 1%引入剩余变量
            a = zeros(m,1);a(j) = -1;
            A = [A a];
            C = [C;0];
    end
end
n1 = size(A,2);%记录转化标准型的未知量的个数
C1 = C;
C = zeros(n1,1);

%找寻单位矩阵
for i = 1:m
    a = 0;
    for j = find(A(i,:)==1)
        if sum(A(:,j)==0) == m-1
            P = [P j];
            a = 1;
        end
    end
    if ~a%若该行无基解,引入人工变量
        j = zeros(m,1);j(i) = 1;
        A = [A j];
        P = [P size(A,2)];
        C = [C;1];
    end
end

P = P(1:m);
CB = C(P);%基变量对应的价值系数
sigma = C'-CB'*inv(A(:,P))*A;
sigma(P) = 0;
SimplexAlgorithmIteration();%第一阶段迭代
if y
    ResultFlag = 0;%如果返回目标值不为0，则不存在最优解
    return;
end
C = C1;
A = A(:,1:size(C,1));
CB = C(P);
sigma = C'-CB'*inv(A(:,P))*A;
sigma(P) = 0;
SimplexAlgorithmIteration();%第二阶段迭代

    function []=SimplexAlgorithmIteration()
        while 1       
            if ~sum(sigma<0)
                if sum(P>n1)%如果基变量含有人工变量
                    return;
                end
                x = zeros(n1,1);
                x(P) = B;
                x = x(1:n);%舍去引入的松弛变量与剩余变量
                if target
                    y = -CB'*B;
                else
                    y = CB'*B;
                end
                ResultFlag = 1;
                return;
            end
            pivot_y = find(sigma==min(sigma));
            pivot_y = pivot_y(1);
            if sum(A(:,pivot_y)<0) == m
                return;
            end
            theta_index = find(A(:,pivot_y)>0);
            theta = B(theta_index)./A(theta_index,pivot_y);
            pivot_x = theta_index(theta==min(theta));%确定主元
            pivot_x = pivot_x(1);
            P(pivot_x) = pivot_y;%更新P
            CB(pivot_x) = C(pivot_y);%更新CB
            %更新系数矩阵
            B(pivot_x) = B(pivot_x)/A(pivot_x,pivot_y);
            A(pivot_x,:) = A(pivot_x,:)./A(pivot_x,pivot_y);
            a = 1:m;
            a(pivot_x) = [];
            for i = a
                B(i) = B(i)-A(i,pivot_y)*B(pivot_x);
                A(i,:) = A(i,:)-A(i,pivot_y)*A(pivot_x,:);    
            end
            sigma = sigma-sigma(pivot_y)*A(pivot_x,:);%更新sigma
        end
    end
end
