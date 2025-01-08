function [R, R_confidence_intervals] = SphereRadius_Sukhovilov3(R0, S, sigma, sigma_m)

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
    
    GAUGE_FREEDOM = 3;
    [X, R, cov_matrix] = optimize_sphere_points2(n, R0, sigma_m, sigma, X0, C0, S, GAUGE_FREEDOM);

    % alpha - significance level (e.g. 0.05 for 95% confidence interval)
    % To obtain confidence intervals at the "3 sigma" level (or 99.73% confidence interval),
    % the significance level alpha should be set so that 99.73% of the distribution
    % is within the interval. This means that alpha will be equal to 1-0.9973=0.0027.
    alpha = 0.0027;
    [CI_lower, CI_upper] = confidence_intervals_ellipsoid([X(:); R], cov_matrix, alpha);
    R_confidence_intervals = [CI_lower(end) CI_upper(end)];
end

function [X, C, R, cov_matrix] = optimize_sphere_points(n, R0, sigma_m, sigma, X0, C0, S, GAUGE_FREEDOM)
    % n - number of points
    % R0 - initial approximation of the sphere radius
    % sigma_m - standard deviation of the error along the radius
    % sigma - standard deviation of the distance measurement error
    % X0 - matrix of initial approximations of point coordinates of size n x 3
    % C0 - row vector of the initial approximation of the coordinates of the center of the sphere

    % Options for the Levenberg-Marquardt algorithm
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'Jacobian', 'on');

    % Initial approximations of the parameters to be optimized
    params0 = [X0(:); C0(:); R0];

    % Function for minimization
    fun = @(params) residuals(params, n, sigma_m, sigma, S);

    % Solve nonlinear least-squares (nonlinear data-fitting) problems
    [params_opt, resnorm, residual, exitflag, output, lambda, J] = lsqnonlin(fun, params0, [], [], options);

    % Extraction of optimal values ??of point coordinates, center and radius
    X = reshape(params_opt(1:n*3), [n, 3]);
    C = params_opt(n*3+1:n*3+3);
    R = params_opt(end);

    % Calculating the covariance matrix
    rank0 = length(params_opt) - GAUGE_FREEDOM;
    cov_matrix = pseudo_inv((J' * J), rank0)*var(residual);
end

function [res, J] = residuals(params, n, sigma_m, sigma,S)
    % Extracting parameters
    X = reshape(params(1:n*3), [n, 3]);
    C = params(n*3+1:n*3+3);
    C = C';
    R = params(end);

    % Calculating distances between points
    D = pdist2(X, X);
    
    % Calculating distances from points to the center of a sphere
    R_dist = sqrt(sum((X - C).^2, 2));

    % Calculating the Jacobian
    J = zeros(n + n*(n-1)/2, n*3 + 3 + 1);
    for i = 1:n
        % Derivatives with respect to coordinates of points
        J(i, (i-1)*3+1:i*3) = (X(i,:) - C) / (R_dist(i) * sigma_m);
        % Derivatives with respect to the coordinates of the center of the sphere
        J(i, n*3+1:n*3+3) = -(X(i,:) - C) / (R_dist(i) * sigma_m);
        % Derivatives with respect to radius
        J(i, end) = -1 / sigma_m;
    end

    % Derivatives with respect to distances between points
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
    
    % Calculating Residuals
    res = [(R_dist - R) / sigma_m; DD / sigma];
   
end

function [X, R, cov_matrix]= optimize_sphere_points2(n, R0, sigma_m, sigma, X0, C0, S, GAUGE_FREEDOM)
    % n - number of points
    % R0 - initial approximation of the sphere radius
    % sigma_m - standard deviation of the error along the radius
    % sigma - standard deviation of the distance measurement error
    % X0 - matrix of initial approximations of point coordinates of size n x 3
    % C0 - row vector of the initial approximation of the coordinates of the center of the sphere
    
    % Translation of coordinates of points into a coordinate system with the center at the center of the sphere
    X0 = X0 - C0';

    % Options for the Levenberg-Marquardt algorithm
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'off', 'Jacobian', 'on');

    % Initial approximations for the parameters to be optimized
    params0 = [X0(:); R0];

    % Function for minimization
    fun = @(params) residuals2(params, n, sigma_m, sigma, S);

    % Solve nonlinear least-squares (nonlinear data-fitting) problems
    [params_opt, resnorm, residual, exitflag, output, lambda, J] = lsqnonlin(fun, params0, [], [], options);

    % Extraction of optimal values ??of point coordinates and radius
    X = reshape(params_opt(1:n*3), [n, 3]);
    R = params_opt(end);
    
    % Calculating the covariance matrix
    rank0 = length(params_opt) - GAUGE_FREEDOM;
    cov_matrix = pseudo_inv((J' * J), rank0)*var(residual);    
end

function [res, J] = residuals2(params, n, sigma_m, sigma, S)
    % Extracting parameters
    X = reshape(params(1:n*3), [n, 3]);
    R = params(end);

    % Calculating distances between points
    D = pdist2(X, X);

    % Calculating distances from points to the center of a sphere
    R_dist = sqrt(sum(X.^2, 2));

    % Calculating the Jacobian
    J = zeros(n + n*(n-1)/2, n*3 + 1);
    for i = 1:n
        % Derivatives with respect to coordinates of points
        J(i, (i-1)*3+1:i*3) = X(i,:) / (R_dist(i) * sigma_m);
        % Derivatives with respect to radius
        J(i, end) = -1 / sigma_m;
    end

    % Derivatives with respect to distances between points
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
        
    % Calculating Residuals
    res = [(R_dist - R) / sigma_m; DD / sigma];
    
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
