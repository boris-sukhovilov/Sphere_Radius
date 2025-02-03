function R = radius_optimal_n_points_on_sphere(S)
    numPoints = size(S,1);
    S2 = 0.5*S.*S;
%     R = sqrt(sum(sum(S2)))/numPoints;
    R = sqrt(kahanSumMatrix(S2))/numPoints;
end
