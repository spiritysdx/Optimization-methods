%这个问题的精确值是x=(1,1)',f(x)=0;
clear all
clc
x0=[-1.2 1]';
[x,val,k]=frcg('fun','gfun',x0);
disp('迭代次数:k=')
disp(k)
disp(['最优解：x = '])
disp(x)
disp(['此时: f(x) = ',num2str(val)]) 

function [x,val,k]=frcg(fun,gfun,x0)
%功能：用共轭梯度法求无约束问题 mini f(x)
%输入：fun,gfun分别是目标函数和梯度，x0是初始点
%输出：x,val分别是近似最优点和最优值，k表示迭代次数
k=0;
maxk=5000;
rho=0.6;
sigma=0.4;
e=1e-5;%精度
n=length(x0);
while(k<maxk)
    g=feval(gfun,x0);%求梯度
    itern=k-(n+1)*floor(k/(n+1));%用于重新开始
     itern=itern+1;
       %计算搜索方向
       if(itern==1)
           d=-g;
       else
           beta=(g'*g)/(g0'*g0);
           d=-g+beta*d0;
           gd=g'*d; %当搜索方向不是下降方向时，插入负梯度方向作为搜索方向
           if(gd>=0.0)
               d=-g;
           end
       end
    if(norm(g)<=e) ,break;end
m=0;
mk=0;
while(m<20)
        if(feval(fun,x0+rho^m*d)<feval(fun,x0)+sigma*rho^m*g'*d);
           mk=m;
           break;
        end
       m=m+1;
end
     x0=x0+d*rho^mk;
    val=feval(fun,x0);
     g0=g;
     d0=d;
     k=k+1;
end
x=x0;
val=feval(fun,x);
end
function f= fun(x)
%目标函数
f=100*(x(1)^2-x(2))^2+(x(1)-1)^2;
end
function  g=gfun(x)
%目标函数的梯度
g=[400*x(1)*(x(1)^2-x(2))+2*(x(1)-1),-200*(x(1)^2-x(2))]';
end


