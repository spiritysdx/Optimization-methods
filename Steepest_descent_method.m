clc
clear
f=@(x1,x2) x1.^2+25*x2.^2;
[x,result]=Min_TD(f,[2;2]);
function [x,result]=Min_TD(f,x0,tol)%fΪ�������������x0Ϊ��ʼ�㣬tolΪ����
if nargin == 2
    tol=1e-6;
end
x01=x0(1);
x02=x0(2);
f_sym=sym(f); %����������ת��Ϊ���ź���
%% F_tdΪ�����ݶ�ֵ  f_tdΪ�˼��㵱ǰ�ݶ�ֵ
F_td=matlabFunction(gradient(f_sym));%����ú����ݶȣ���ȡ���ݶȺ������
f_td=F_td(x01,x02); %���㵱ǰ���ݶ�ֵ
d_k=-f_td/norm(f_td);%�½�����
%% Ѱ����ѵ������� alfa
while norm(f_td) > tol %��ǰ�ݶ�ֵ������Ҫ��
    syms alfa
    x1=x0(:)+alfa*d_k;%����������������ʽ����alfa 
    x11=x1(1);
    x12=x1(2);
    fx1=f(x11,x12);%����x1������ԭ�������ʽ
    d_x1=diff(fx1);%��ԭ���������ʽ��
    d_alfa=double(solve(d_x1));%�����ʽΪ0ʱ��alfa������Ϊ��ʱ�ĵ�Ϊ��Сֵ�㣩
    
    x0=x0(:)+d_alfa*d_k;%���е���
    x01=x0(1);
    x02=x0(2);
    f_td=F_td(x01,x02);%��ǰ�ݶ�ֵ
    if norm(f_td) < tol %����Ҫ���˳�
        break;
    end
    d_k=-f_td/norm(f_td);%���򣬼�������
end
x=x0;
result=f(x(1),x(2));
end

 
