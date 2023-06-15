clc
clear
f=@(x1,x2) x1.^2+25*x2.^2;
[x,result]=Min_TD(f,[2;2]);
function [x,result]=Min_TD(f,x0,tol)%f为匿名函数句柄，x0为初始点，tol为精度
if nargin == 2
    tol=1e-6;
end
x01=x0(1);
x02=x0(2);
f_sym=sym(f); %将匿名函数转化为符号函数
%% F_td为计算梯度值  f_td为了计算当前梯度值
F_td=matlabFunction(gradient(f_sym));%计算该函数梯度，并取得梯度函数句柄
f_td=F_td(x01,x02); %计算当前的梯度值
d_k=-f_td/norm(f_td);%下降方向
%% 寻找最佳迭代步长 alfa
while norm(f_td) > tol %当前梯度值不满足要求
    syms alfa
    x1=x0(:)+alfa*d_k;%将这个点带入迭代表达式，求alfa 
    x11=x1(1);
    x12=x1(2);
    fx1=f(x11,x12);%计算x1带入后的原函数表达式
    d_x1=diff(fx1);%对原函数表达是式求导
    d_alfa=double(solve(d_x1));%求解表达式为0时的alfa（导数为零时的点为极小值点）
    
    x0=x0(:)+d_alfa*d_k;%进行迭代
    x01=x0(1);
    x02=x0(2);
    f_td=F_td(x01,x02);%当前梯度值
    if norm(f_td) < tol %满足要求，退出
        break;
    end
    d_k=-f_td/norm(f_td);%否则，继续迭代
end
x=x0;
result=f(x(1),x(2));
end

 
