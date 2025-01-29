% ------------------------------------------------------------------------------------------
% Calculate the coordinates of the vectors of n points in the coordinate system
% with the center at the center of gravity of the coordinates of the vectors of points
% and the axes directed along the eigenvectors, corresponding to the 3 largest eigenvalues
% matrix of scalar products of the coordinates of the vectors of points
% Input:
% s - matrix of pairwise distances
% Output:
% cg - matrix of vector coordinates of size 3xn
% ------------------------------------------------------------------------------------------
function [cg] = CenterOfGravityCoordFromPairDistance(s)
    n = size(s,1);

    % Calculation of the matrix of scalar products in the coordinate system
    % with the center at the center of gravity
    J=eye(n,n)-(1/n)*ones(n,n);
    c = -0.5*J*s.^2*J;

    % matrix of vector coordinates of size 3xn
    [cg]=factorization2(c);
end

% ---------------------------------------------------------------------------------------
% calculation of coordinates, decomposition of the matrix of scalar products via svd
% ---------------------------------------------------------------------------------------
function [beta]=factorization2(c)
    % For a square symmetric matrix, the singular values are equal to the absolute values of the eigenvalues
    [U,L,~] = svd(c);
    di = diag(L);

    % we take 3 vectors corresponding to the 3 largest singular values
    U = U(:,1:3);
    sqrt_lamda = diag(sqrt(di(1:3)));

    beta = sqrt_lamda*U';
end
