%内点法
syms x1 x2
f=x1^3+x2^2;
h=x1+x2-1;
[X1,X2]=IPM(f,h);
function [X1, X2] = IPM(f, h)
%       IPM     内点法
%       f       目标函数
%       h       等式约束
R = 1;          %障碍因子R值
C = 0.1;        %倍率
eps =  0.001;   %误差
syms x1 x2;             %定义变量

H = matlabFunction(log(h));%函数句柄
F = matlabFunction(f);%函数句柄
while true
    G = f - R * log(h);    %构造函数
    G_x1 = diff(G,x1);      %求x1偏导
    G_x2 = diff(G,x2);      %求x2偏导
    X = solve(G_x1==0,G_x2==0,x1,x2); %解方程组
    X1 = real(double(X.x1));
    X2 = real(double(X.x2));
    if abs(R * H(X1, X2)) <= eps     % 判断终止条件
        break
    end
    R = R * C;
end
%%%%%%%%%%%% 多个最优解则选择目标函数最小的
i = 1;
answer = 10000000;
X = [X1 X2];
[n, ~] = size(X);
for j = 1:1:n
    if F(X(j,1),X(j,2)) < answer && X(j,1)>0 && X(j,2)>0
        answer = F(X(j,1),X(j,2));
        i = j;
    end
end
%%%%%%%%%%% 更新答案
X1 = X1(i);
X2 = X2(i);
end
