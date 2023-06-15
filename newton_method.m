%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%文件名:Newton.m
%
%f表示函数表达式
%H0表示初始的海森矩阵
%x0表示初始的迭代点 为列向量
%m表示变量的个数
%k表示迭代次数
%X存储每次迭代的x,F为函数值，G为每次的梯度，H为海森阵，HN为海森矩阵的逆
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x1; % 变量x1
syms x2; % 变量x2
f = (x1 - 3).^4 + (x1 - 3*x2).^2; % 函数表达式
x0=[0 0]'; % 初始迭代点
H0=[110 -6;-6 18]; % 初始海森矩阵
m=2; % 变量个数
k=30; % 迭代次数

[X, H, F, G, HN] = Newton(f,H0,x0,m,k);

function[X, H, F, G,HN] = Newton(f,H0,x0,m,k)
    x1 = sym('x',[1,m]); % [x1, x2]
    %f = (x1(1) - 3)^4 + (x1(1) - 3*x1(2))^2;
    c = num2cell(x1); % c=变量[x1, x2]
    g = sym('x',[m,1]); % [x1, x2]^T
    
    X = zeros(m, k+1); % x1、x2的迭代值
    H = zeros( m, m, k+1); % hessian的迭代值
    F = zeros(1, k+1); % function的迭代值
    G = zeros(m, k+1); % function‘的迭代值
    HN = zeros( m, m, k+1); % hessian的逆阵的迭代值
    
    H(:,:,1) = H0; % hessian初始化
    HN(:,:,1) = inv(H0); % hessian逆初始化
    X(:,1) = x0; % X(x1, x2)初始化
    F(1,1) = subs(f, c, {X(:,1)'}); % 初始X值赋予F
    h = hessian(f,x1);%求海森矩阵
    
    for n = 1:m % f对x1、x2分别求偏导
        g(n) = diff(f,x1(n));
    end
     G(:,1) = subs(g,c,{X(:,1)'}); % 初始X导赋予G
     
     % 迭代
    for n = 1:k
        X(:,n+1) = X(:,n) - (H(:,:,n))\G(:,n);
        F(1,n+1) = subs(f,c,{X(:,n+1)'}); 
        G(:,n+1) = subs(g,c,{X(:,n+1)'});
        H(:,:,n+1) = subs(h,c,{X(:,n+1)'});
        HN(:,:,n+1) = inv(H(:,:,n+1));
    end
end
