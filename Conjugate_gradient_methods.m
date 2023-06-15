%�������ľ�ȷֵ��x=(1,1)',f(x)=0;
clear all
clc
x0=[-1.2 1]';
[x,val,k]=frcg('fun','gfun',x0);
disp('��������:k=')
disp(k)
disp(['���Ž⣺x = '])
disp(x)
disp(['��ʱ: f(x) = ',num2str(val)]) 

function [x,val,k]=frcg(fun,gfun,x0)
%���ܣ��ù����ݶȷ�����Լ������ mini f(x)
%���룺fun,gfun�ֱ���Ŀ�꺯�����ݶȣ�x0�ǳ�ʼ��
%�����x,val�ֱ��ǽ������ŵ������ֵ��k��ʾ��������
k=0;
maxk=5000;
rho=0.6;
sigma=0.4;
e=1e-5;%����
n=length(x0);
while(k<maxk)
    g=feval(gfun,x0);%���ݶ�
    itern=k-(n+1)*floor(k/(n+1));%�������¿�ʼ
     itern=itern+1;
       %������������
       if(itern==1)
           d=-g;
       else
           beta=(g'*g)/(g0'*g0);
           d=-g+beta*d0;
           gd=g'*d; %�������������½�����ʱ�����븺�ݶȷ�����Ϊ��������
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
%Ŀ�꺯��
f=100*(x(1)^2-x(2))^2+(x(1)-1)^2;
end
function  g=gfun(x)
%Ŀ�꺯�����ݶ�
g=[400*x(1)*(x(1)^2-x(2))+2*(x(1)-1),-200*(x(1)^2-x(2))]';
end


