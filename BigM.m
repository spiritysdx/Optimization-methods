function [ x,y ] = BigM( f,A,b )
%����f�Ǽ����������飬1*nά
%����A��Լ������ m*nά
%����b��Լ�������� 1*mά
%���x�ǽ�����
%���y�����Ž�
%�ж�����ά���Ƿ����
%����ʼ�����α�,����M����
[n,m]=size(A);%n��m��
M=10000;
S=[f -1*M*ones(1,n) 0;
   A   eye(n)       b'];
format rat %������Է�����ʾ
[n,m]=size(S);
%���˹������ļ���������
for k=1:n-1
    S(1,:)=S(1,:)+S(k+1,:)*M;   
end
%�жϼ����� r<=0
r=find(S(1,1:m-1)>0);
len=length(r);
flag=0;
%�д���0�ļ����������ѭ��
while(len)
    %���Ǹ�����������������Ԫ���Ƿ�С�ڵ���0
    for k=1:length(r)
        d=find(S(:,r(k))>0);
        if(length(d)+1==2)
        error('�����Ž⣡����')  
        %break;
        end
    end
    %�ҵ������������ֵ
    [Rk,j]=max(S(1,1:m-1));
    %���ֵ�����б�ֵΪ��������Сֵbr/a_rk
    br=S(2:n,m)./S(2:n,j);
    %�ѱ�ֵ�еĸ�����������
    for p=1:length(br)
        if(br(p)<0)br(p)=Inf;
        end
    end
    [h,i]=min(br);%��������ֵ��Сֵ
    % i+1Ϊת��Ԫ�к�(��S��)��jΪת��Ԫ�к�
    i=i+1;
    %���л���,ת��Ԫ��1
    S(i,:)=S(i,:)./S(i,j);
    %ת��Ԫ����������Ԫ�ض���0
    for k=1:n
        if(k~=i)
            S(k,:)=S(k,:)-S(i,:)*S(k,j);
        end   
    end
    %�жϼ����� r<=0
    r=find(S(1,1:m-1)>0);
    len=length(r);
%     %�����ã�����ѭ������
%     if(len>0)flag=flag+1;end
%     if(flag==2)break;end
%     S
end
    %������ȫ���������ҵ����Ž�
    %�ǻ�������0
    x=zeros(1,m-1);
    for i=1:m-1
        %�ҵ�������
        j=find(S(:,i)==1);%ÿ����1�ĸ���
        k=find(S(:,i)==0);%ÿ����0�ĸ���
        if((length(j)+1==2)&&(length(k)+1==n))
            %iΪ����Ԫ�кţ�j���к�
            x(i)=S(j,m);
        end
    end
    y=S(1,m);%���Ž�
    S
end