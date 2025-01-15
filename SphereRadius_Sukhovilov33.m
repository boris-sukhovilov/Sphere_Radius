function [R, R_confidence_intervals] = SphereRadius_Sukhovilov33(R0, S, sigma, sigma_m)

    n = size(S,1);
    
    % Initial approximations for point coordinates
    [cg] = CenterOfGravityCoordFromPairDistance(S);
    X0 = cg';
    C0 = compute_center(X0);
    
%     r = 0;
%     for i = 1 : n
%         r = r + norm(cg(:,i)-C0);
%     end
%     r = r/n;    
%     r
    
%     GAUGE_FREEDOM = 3;
%     [X, R, cov_matrix] = optimize_sphere_points2(n, R0, sigma_m, sigma, X0, C0, S, GAUGE_FREEDOM);
    
    [R, theta, phi] = optimize_spherical_coords(X0, S, sigma, sigma_m, R0, C0);
    R
    R_confidence_intervals = [];
    
    % alpha - significance level (e.g. 0.05 for 95% confidence interval)
    % To obtain confidence intervals at the "3 sigma" level (or 99.73% confidence interval),
    % the significance level alpha should be set so that 99.73% of the distribution
    % is within the interval. This means that alpha will be equal to 1-0.9973=0.0027.
%     alpha = 0.0027;
%     [CI_lower, CI_upper] = confidence_intervals_ellipsoid([X(:); R], cov_matrix, alpha);
%     R_confidence_intervals = [CI_lower(end) CI_upper(end)];
end

function [R, theta, phi] = optimize_spherical_coords(points, d, sigma, sigma_m, R0, C0)
    % Начальные приближения декартовых координат
    x = points(:, 1) - C0(1);
    y = points(:, 2) - C0(2);
    z = points(:, 3) - C0(3);

    % Начальные приближения сферических координат и радиуса
    n = size(points, 1);
    R = R0;
    theta = acos(z ./ sqrt(x.^2 + y.^2 + z.^2));
%     theta = acos(z ./ R);
    phi = atan2(y, x);

    % Опции для оптимизации
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'Jacobian', 'on','TolX',1.0e-15,'TolFun',1.0e-15);

    % Функция для минимизации
    fun = @(params) residuals(params, n, d, sigma, sigma_m);

    % Начальные параметры
    params0 = [R; theta; phi];

    % Минимизация функционала
    params_opt = lsqnonlin(fun, params0, [], [], options);

    % Разделение оптимизированных параметров
    R = params_opt(1);
    theta = params_opt(2:n+1);
    phi = params_opt(n+2:2*n+1);
end

function [res, J] = residuals(params, n, d, sigma, sigma_m)
    R = params(1);
    theta = params(2:n+1);
    phi = params(n+2:2*n+1);

    % Вычисление декартовых координат из сферических
    x = R .* sin(theta) .* cos(phi);
    y = R .* sin(theta) .* sin(phi);
    z = R .* cos(theta);

    % Инициализация остатков и якобиана
    res = zeros(n*(n-1)/2, 1);
    J = ones(n*(n-1)/2, 2*n+1);

    k = 1;
    for i = 1:n-1
        for j = i+1:n
            d_ij = sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2 + (z(i) - z(j))^2);
            res(k) = (d(i,j) - d_ij) / sigma;

            % Вычисление частных производных для якобиана
            dx_dR_i = sin(theta(i)) * cos(phi(i));
            dy_dR_i = sin(theta(i)) * sin(phi(i));
            dz_dR_i = cos(theta(i));
            dx_dtheta_i = R * cos(theta(i)) * cos(phi(i));
            dy_dtheta_i = R * cos(theta(i)) * sin(phi(i));
            dz_dtheta_i = -R * sin(theta(i));
            dx_dphi_i = -R * sin(theta(i)) * sin(phi(i));
            dy_dphi_i = R * sin(theta(i)) * cos(phi(i));
            dz_dphi_i = 0;

            dx_dtheta_j = R * cos(theta(j)) * cos(phi(j));
            dy_dtheta_j = R * cos(theta(j)) * sin(phi(j));
            dz_dtheta_j = -R * sin(theta(j));
            dx_dphi_j = -R * sin(theta(j)) * sin(phi(j));
            dy_dphi_j = R * sin(theta(j)) * cos(phi(j));
            dz_dphi_j = 0;

            J(k, 1) = ((x(i) - x(j)) * dx_dR_i + (y(i) - y(j)) * dy_dR_i + (z(i) - z(j)) * dz_dR_i) / (d_ij * sigma);
            J(k, i+1) = ((x(i) - x(j)) * dx_dtheta_i + (y(i) - y(j)) * dy_dtheta_i + (z(i) - z(j)) * dz_dtheta_i) / (d_ij * sigma);
            J(k, n+i+1) = ((x(i) - x(j)) * dx_dphi_i + (y(i) - y(j)) * dy_dphi_i + (z(i) - z(j)) * dz_dphi_i) / (d_ij * sigma);
            J(k, j+1) = -((x(i) - x(j)) * dx_dtheta_j + (y(i) - y(j)) * dy_dtheta_j + (z(i) - z(j)) * dz_dtheta_j) / (d_ij * sigma);
            J(k, n+j+1) = -((x(i) - x(j)) * dx_dphi_j + (y(i) - y(j)) * dy_dphi_j + (z(i) - z(j)) * dz_dphi_j) / (d_ij * sigma);

            k = k + 1;
        end
    end

%     % Добавление остатков и якобиана для погрешностей вдоль радиуса сферы
%     for i = 1:n
%         r_i = sqrt(x(i)^2 + y(i)^2 + z(i)^2);
%         res(k) = (R - r_i) / sigma_m;
% 
% %         dx_dR_i = sin(theta(i)) * cos(phi(i));
% %         dy_dR_i = sin(theta(i)) * sin(phi(i));
% %         dz_dR_i = cos(theta(i));
% 
%         dx_dtheta_i = R * cos(theta(i)) * cos(phi(i));
%         dy_dtheta_i = R * cos(theta(i)) * sin(phi(i));
%         dz_dtheta_i = -R * sin(theta(i));
%         dx_dphi_i = -R * sin(theta(i)) * sin(phi(i));
%         dy_dphi_i = R * sin(theta(i)) * cos(phi(i));
%         dz_dphi_i = 0;
%             
%         J(k, 1) = (R - r_i) / (r_i * sigma_m);
%         J(k, i+1) = (x(i) * dx_dtheta_i + y(i) * dy_dtheta_i + z(i) * dz_dtheta_i) / (r_i * sigma_m);
%         J(k, n+i+1) = (x(i) * dx_dphi_i + y(i) * dy_dphi_i + z(i) * dz_dphi_i) / (r_i * sigma_m);
% 
%         k = k + 1;
%     end
    
end


function C = compute_center(X)
    % X is a matrix of point coordinates of size n x 3

    % Number of points
    n = size(X, 1);

    A = zeros(n-1, 3);
    b = zeros(n-1, 1);

    for i = 2:n
        A(i-1, :) = 2 * (X(i, :) - X(1, :));
        b(i-1) = sum(X(i, :).^2 - X(1, :).^2);
    end

    C = A \ b;
end
