function [ x,y ] = BigM( f,A,b )
%输入f是检验数的数组，1*n维
%输入A是约束矩阵， m*n维
%输入b是约束向量， 1*m维
%输出x是解向量
%输出y是最优解
%判断输入维数是否相符
%做初始单纯形表,加入M变量
[n,m]=size(A);%n行m列
M=10000;
S=[f -1*M*ones(1,n) 0;
   A   eye(n)       b'];
format rat %将结果以分数表示
[n,m]=size(S);
%将人工变量的检验数置零
for k=1:n-1
    S(1,:)=S(1,:)+S(k+1,:)*M;   
end
%判断检验数 r<=0
r=find(S(1,1:m-1)>0);
len=length(r);
flag=0;
%有大于0的检验数则进入循环
while(len)
    %检查非负检验数所对列向量元素是否都小于等于0
    for k=1:length(r)
        d=find(S(:,r(k))>0);
        if(length(d)+1==2)
        error('无最优解！！！')  
        %break;
        end
    end
    %找到检验数中最大值
    [Rk,j]=max(S(1,1:m-1));
    %最大值所在列比值为正数且最小值br/a_rk
    br=S(2:n,m)./S(2:n,j);
    %把比值中的负数都变无穷
    for p=1:length(br)
        if(br(p)<0)br(p)=Inf;
        end
    end
    [h,i]=min(br);%列向量比值最小值
    % i+1为转轴元行号(在S中)，j为转轴元列号
    i=i+1;
    %进行换基,转轴元置1
    S(i,:)=S(i,:)./S(i,j);
    %转轴元所在列其他元素都置0
    for k=1:n
        if(k~=i)
            S(k,:)=S(k,:)-S(i,:)*S(k,j);
        end   
    end
    %判断检验数 r<=0
    r=find(S(1,1:m-1)>0);
    len=length(r);
%     %调试用，控制循环步数
%     if(len>0)flag=flag+1;end
%     if(flag==2)break;end
%     S
end
    %检验数全部非正，找到最优解
    %非基变量置0
    x=zeros(1,m-1);
    for i=1:m-1
        %找到基变量
        j=find(S(:,i)==1);%每列中1的个数
        k=find(S(:,i)==0);%每列中0的个数
        if((length(j)+1==2)&&(length(k)+1==n))
            %i为基本元列号，j是行号
            x(i)=S(j,m);
        end
    end
    y=S(1,m);%最优解
    S
end