function [R, sigma_R] = SphereRadius_Sukhovilov3(R0, S, sigma, sigma_m)

    n = size(S,1);
    
    % ��������� ����������� ��� ��������� ����� � ���� ���� ����������
    [cg] = CenterOfGravityCoordFromPairDistance(S);
    X0 = cg';
    C0 = compute_center(X0, R0);
    
    [X, C, R, cov_matrix] = optimize_sphere_points(n, R0, sigma_m, sigma, X0, C0, S);
    X, C, R
%     cov_matrix
    
    svd(cov_matrix)
   
    sigma_R = sqrt(cov_matrix(end,end));
end

function [X, C, R, cov_matrix] = optimize_sphere_points(n, R0, sigma_m, sigma, X0, C0, S)
    % n - ���������� �����
    % R0 - ��������� ����������� ������� �����
    % sigma_m - ��� ����������� ����� �������
    % sigma - ��� ����������� ��������� ����������
    % X0 - ��������� ����������� ��������� ����� (n x 3)
    % C0 - ��������� ����������� ��������� ������ ����� (1 x 3)

    % ����� ��� ��������� ����������-����������
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'Jacobian', 'on');

    % ��������� ����������� ��� �������������� ����������
    params0 = [X0(:); C0(:); R0];

    % ������� ��� �����������
    fun = @(params) residuals(params, n, sigma_m, sigma, S);

    % ����������� �����������
    [params_opt, resnorm, residual, exitflag, output, lambda, J] = lsqnonlin(fun, params0, [], [], options);

    % ���������� ����������� �������� ��������� �����, ������ � �������
    X = reshape(params_opt(1:n*3), [n, 3]);
    C = params_opt(n*3+1:n*3+3);
    R = params_opt(end);

    % ���������� ������� ����������
%     cov_matrix = inv(J' * J) * var(residual);
    cov_matrix = (J' * J)\eye(length(params_opt))*var(residual);
end

function [res, J] = residuals(params, n, sigma_m, sigma,S)
    % ���������� ����������
    X = reshape(params(1:n*3), [n, 3]);
    C = params(n*3+1:n*3+3);
    C = C';
    R = params(end);

    % ���������� ���������� ����� �������
    D = pdist2(X, X);
    
    % ���������� ���������� �� ����� �� ������ �����
    R_dist = sqrt(sum((X - C).^2, 2));

    % ���������� ��������
    J = zeros(n + n*(n-1)/2, n*3 + 3 + 1);
    for i = 1:n
        % ����������� �� ����������� �����
        J(i, (i-1)*3+1:i*3) = (X(i,:) - C) / (R_dist(i) * sigma_m);
        % ����������� �� ����������� ������
        J(i, n*3+1:n*3+3) = -(X(i,:) - C) / (R_dist(i) * sigma_m);
        % ����������� �� �������
        J(i, end) = -1 / sigma_m;
    end

    % ����������� �� ����������� ����� �������
    DD = zeros(n*(n-1)/2,1);
    k = 1;
    idx = n;
    for i = 1:n-1
        for j = i+1:n
            idx = idx + 1;
            J(idx, (i-1)*3+1:i*3) = (X(i,:) - X(j,:)) / (D(i,j) * sigma);
            J(idx, (j-1)*3+1:j*3) = -(X(i,:) - X(j,:)) / (D(i,j) * sigma);
                        
            DD(k) = D(i,j) - S(i,j);
        end
    end
    
    % ���������� ��������
    res = [(R_dist - R) / sigma_m; DD / sigma];
   
end

function C = compute_center(X, R)
    % X - ������� ��������� ����� (n x 3)
    % R - ������ �����

    % ���������� �����
    n = size(X, 1);

    % ����������� ������� ���������
    A = zeros(n-1, 3);
    b = zeros(n-1, 1);

    for i = 2:n
        A(i-1, :) = 2 * (X(i, :) - X(1, :));
        b(i-1) = sum(X(i, :).^2 - X(1, :).^2);
    end

    % ������� ������� �������� ���������
    C = A \ b;
end
