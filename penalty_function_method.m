%�ڵ㷨
syms x1 x2
f=x1^3+x2^2;
h=x1+x2-1;
[X1,X2]=IPM(f,h);
function [X1, X2] = IPM(f, h)
%       IPM     �ڵ㷨
%       f       Ŀ�꺯��
%       h       ��ʽԼ��
R = 1;          %�ϰ�����Rֵ
C = 0.1;        %����
eps =  0.001;   %���
syms x1 x2;             %�������

H = matlabFunction(log(h));%�������
F = matlabFunction(f);%�������
while true
    G = f - R * log(h);    %���캯��
    G_x1 = diff(G,x1);      %��x1ƫ��
    G_x2 = diff(G,x2);      %��x2ƫ��
    X = solve(G_x1==0,G_x2==0,x1,x2); %�ⷽ����
    X1 = real(double(X.x1));
    X2 = real(double(X.x2));
    if abs(R * H(X1, X2)) <= eps     % �ж���ֹ����
        break
    end
    R = R * C;
end
%%%%%%%%%%%% ������Ž���ѡ��Ŀ�꺯����С��
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
%%%%%%%%%%% ���´�
X1 = X1(i);
X2 = X2(i);
end
