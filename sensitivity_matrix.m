% Calculation of the sensitivity matrix for the solution R from the implicit equation: f(R) = det(R^2*ones(n)-S.*S/2) = 0
% sensitivity_matrix = - (df/dS)/(df/dR)
% Automatic differentiation tools in matrix form: https://www.matrixcalculus.org/matrixCalculus
% df/dS = -(M').*S, где M = adjointMatrix(R^2*ones(n)-S.*S/2);
% df/dR =  2*R*trace(ones(n)*M);
function sm = sensitivity_matrix(R, S)
    n = size(S,1);
    M = adjointMatrix(R^2*ones(n)-S.*S/2);
    df_dS = -(M').*S;
    df_dR =  2*R*trace(ones(n)*M);
    sm = -df_dS./df_dR;
end