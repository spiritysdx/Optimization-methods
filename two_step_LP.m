A=[1 -2 1;4 -1 -2;-2 0 1];
B=[11 -3 1];
B=B';
C=[3 -1 -1];
C=C';
target=1;
sign=[-1 -1 0];
sign=sign';
[x,y,result]=TwoStageSimplexAlgorithm(A,B,C,target,sign);
function [x,y,ResultFlag]=TwoStageSimplexAlgorithm(A,B,C,target,sign)
% 2020-4-2 ��orz
%inputs:
%   A:ϵ������ m*n
%   B:�Ҷ����� m*1
%   C:�۸�ϵ������ n*1
%alternative inputs:
%   target:�Ż�Ŀ�� 0 ~ min; 1 ~ max;
%   sign:Լ�������������� m*1 ����
%       -1 ~ '<=';0 ~ '=';1 ~ '>'; Ĭ��Ϊ0,��Ϊ��ʽ
%outputs:
%   x:���Ž� n*1
%   y:����ֵ num
%   ResultFlag:�Ƿ��ҵ����Ž�
[m,n] = size(A);
P = [];
x = zeros(n,1);
y = 0;
ResultFlag = 0;
j = 0;

%standardization
if target
    C = -C;%Ŀ�꺯����ת��
end
A(B<0,:) = -A(B<0,:);
sign(B<0,:) = -sign(B<0,:);
B = abs(B);%Լ��������ת��
for i = sign'
    j = j+1;
    switch i
        case -1%�����ɳڱ���
            a = zeros(m,1);a(j) = 1;
            A = [A a];
            C = [C;0];
        case 1%����ʣ�����
            a = zeros(m,1);a(j) = -1;
            A = [A a];
            C = [C;0];
    end
end
n1 = size(A,2);%��¼ת����׼�͵�δ֪���ĸ���
C1 = C;
C = zeros(n1,1);

%��Ѱ��λ����
for i = 1:m
    a = 0;
    for j = find(A(i,:)==1)
        if sum(A(:,j)==0) == m-1
            P = [P j];
            a = 1;
        end
    end
    if ~a%�������޻���,�����˹�����
        j = zeros(m,1);j(i) = 1;
        A = [A j];
        P = [P size(A,2)];
        C = [C;1];
    end
end

P = P(1:m);
CB = C(P);%��������Ӧ�ļ�ֵϵ��
sigma = C'-CB'*inv(A(:,P))*A;
sigma(P) = 0;
SimplexAlgorithmIteration();%��һ�׶ε���
if y
    ResultFlag = 0;%�������Ŀ��ֵ��Ϊ0���򲻴������Ž�
    return;
end
C = C1;
A = A(:,1:size(C,1));
CB = C(P);
sigma = C'-CB'*inv(A(:,P))*A;
sigma(P) = 0;
SimplexAlgorithmIteration();%�ڶ��׶ε���

    function []=SimplexAlgorithmIteration()
        while 1       
            if ~sum(sigma<0)
                if sum(P>n1)%��������������˹�����
                    return;
                end
                x = zeros(n1,1);
                x(P) = B;
                x = x(1:n);%��ȥ������ɳڱ�����ʣ�����
                if target
                    y = -CB'*B;
                else
                    y = CB'*B;
                end
                ResultFlag = 1;
                return;
            end
            pivot_y = find(sigma==min(sigma));
            pivot_y = pivot_y(1);
            if sum(A(:,pivot_y)<0) == m
                return;
            end
            theta_index = find(A(:,pivot_y)>0);
            theta = B(theta_index)./A(theta_index,pivot_y);
            pivot_x = theta_index(theta==min(theta));%ȷ����Ԫ
            pivot_x = pivot_x(1);
            P(pivot_x) = pivot_y;%����P
            CB(pivot_x) = C(pivot_y);%����CB
            %����ϵ������
            B(pivot_x) = B(pivot_x)/A(pivot_x,pivot_y);
            A(pivot_x,:) = A(pivot_x,:)./A(pivot_x,pivot_y);
            a = 1:m;
            a(pivot_x) = [];
            for i = a
                B(i) = B(i)-A(i,pivot_y)*B(pivot_x);
                A(i,:) = A(i,:)-A(i,pivot_y)*A(pivot_x,:);    
            end
            sigma = sigma-sigma(pivot_y)*A(pivot_x,:);%����sigma
        end
    end
end
