
clc
clear
f=@(x) x.^2-x+2;  %����һ����������
f1=@(x) 2*x-1;   %����һ������������Ϊf�ĵ�����
x1=vpa(solve(f1),5); %��ֵ��
[k,x,result]=fibonacci(f,[-1 3],1e-9)
tol=x-double(x1(1))   %�����

function [k,x,result]=fibonacci(f,interval,delta)   %kΪ��������
a=interval(1);
b=interval(2);
F=[1 1];      %쳲������������еĵ�һ��͵ڶ���
n=3;
while F(end)<(b-a)/delta    %����쳲���������
    F(n)=F(n-1)+F(n-2);    
    n=n+1;
end
n=n-1;  %����쳲��������е�ʱ��ĩβ�����һ��1���˴�Ҫ��ȥһ��1
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
x=(a+b)/2;     %ȡ����������е�Ϊ��Сֵ��
result=f(x);
end
