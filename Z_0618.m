%%%%      0.618��      %%%%%
%%%%%    Liu Deping    %%%%%
%%%%%    2020.09.21    %%%%%
clc;
clear all
format short
s=input('�����뺯�����ʽ��f = ','s');
f=inline(s);
a = input('������������˵�a��ֵ:');
b = input('������������˵�b��ֵ:');
eps= input('���������㾫��:');    %%��|b-a|<=eps��
k=0;
x1=a+0.382*(b-a);
x2=a+0.618*(b-a);
fprintf(' k          [a,b]          x1        x2        f(x1)       f(x2)\n ');
fprintf('%d      [%.3f,%.3f]    %.3f     %.3f     %.3f      %.3f\n', k,a,b,x1,x2,f(x1),f(x2));
while abs(a-b)>eps
    k=k+1;
    x1=a+0.382*(b-a);
    x2=a+0.618*(b-a);
    if  f(x1)>f(x2)
        a=x1;
        b=b;
    elseif  f(x1)<=f(x2)
        a=a;
        b=x2;
    end
   fprintf(' %d      [%.3f,%.3f]    %.3f     %.3f     %.3f      %.3f\n', k,a,b,x1,x2,f(x1),f(x2));
end
fprintf('����%d�ε�������ȡ%.3f��Ϊ�������Ž�',k,(a+b)/2);

