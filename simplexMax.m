function [x, result_case] = simplexMax(c, A, b, x0)
 
ind_B = find(x0 ~=0);
ind_N = find(x0 == 0);
 
while (1)
    B = A(:, ind_B);
    N = A(:, ind_N);
    x_B = x0(ind_B);
    x_N = x0(ind_N);
    c_B = c(ind_B);
    c_N = c(ind_N);
 
    s = c_N - N' * inv(B)' * c_B; 
    %'
    if isempty(find(sign(s) > 0, 1))
        % found optimal solution
        x = x0;
        result_case = 0;
        return
    end
    % ��֤ s ������
    % ѡ��������������Ϊ�������� q
    [max_q i_q] = max(s);
    q = ind_N(i_q);
    d = B \ A(:, q);
    if isempty(find(sign(d) > 0, 1))
        % unbound case
        x = [];
        result_case = 1;
        break;
    end
    % ��֤ d ������
    % ѡ����С�Ǹ�ratio ��Ϊ������� p
    ratio_array = x_B./d;
    ratio = min(ratio_array((ratio_array > 0)));
 
    % ���� x0 ��ֵ
    x0(ind_B) = x_B - ratio * d;
    e_iq = zeros(length(x_N), 1);
    e_iq(i_q) = 1;
    x0(ind_N) = x_N + ratio * e_iq;
    % ���»��������ǻ���������
    ind_B = find(x0 ~=0);
    ind_N = find(x0 == 0);
end
end