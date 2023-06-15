
clc
clear
f=@(x) x.^2-x+2;  %定义一个匿名函数
f1=@(x) 2*x-1;   %定义一个匿名函数，为f的导函数
x1=vpa(solve(f1),5); %求极值点
[k,x,result]=fibonacci(f,[-1 3],1e-9)
tol=x-double(x1(1))   %求误差

function [k,x,result]=fibonacci(f,interval,delta)   %k为迭代次数
a=interval(1);
b=interval(2);
F=[1 1];      %斐波那契数列数列的第一项和第二项
n=3;
while F(end)<(b-a)/delta    %产生斐波那契数列
    F(n)=F(n-1)+F(n-2);    
    n=n+1;
end
n=n-1;  %产生斐波那契数列的时候末尾多加了一个1，此处要减去一个1
x1=a+(F(n-2)/F(n))*(b-a);
x2=a+(F(n-1)/F(n))*(b-a);
f1=f(x1);
f2=f(x2);
k=0;
while n>3 
     k=k+1;
     n=n-1;
    if f1>f2
        a=x1;
        x1=x2;
        f1=f2;
        x2=a+(F(n-1)/F(n))*(b-a);
        f2=f(x2);
    elseif f1==f2
        a=x1;
        b=x2;
        x1=a+(F(n-2)/F(n))*(b-a);
        x2=a+(F(n-1)/F(n))*(b-a);
    else
        b=x2;
        x2=x1;
        f2=f1;
        x1=a+(F(n-2)/F(n))*(b-a);
        f1=f(x1);
    end
end
x=(a+b)/2;     %取最终区间的中点为极小值点
result=f(x);
end
