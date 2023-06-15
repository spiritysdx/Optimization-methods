%% 第一阶段
A1 = [1 1 -1 0 0 1 0; 1 0 0 -1 0 0 1; 2 1 0 0 1 0 0];
b1 = [350; 125; 600];
c1 = [0 0 0 0 0 -1 -1];
ind1 = [6 7 5];
[x1, z1, ST1, ca1] = simplexMax(c1, A1, b1, ind1)

%% 第二阶段
A2 = ST1(end-m1:end-1,4:end-3);   %2个人工变量, 2+1
b2 = ST1(end-m1:end-1,3);
c2 = [-2 -3 0 0 0];
ind2 = ST1(end-m1:end-1,2)';
[x2, z2, ST2, ca2] = simplexMax(c2, A2, b2, ind2)

%% https://www.codeleading.com/article/87071107518/