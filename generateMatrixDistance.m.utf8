function [S, S0] = generateMatrixDistance(points, sigma)

    numPoints = size(points, 2);
    
    S0 = zeros(numPoints,numPoints);
    S = S0;
    
    for i = 1 : numPoints - 1
        for j = i + 1 : numPoints
            S0(i,j) = norm(points(:,i) - points(:,j));
            S0(j,i) = S0(i,j);
            err=sigma*randn;
            S(i,j)= S0(i,j)+err;
            S(j,i)=S(i,j);
        end
    end
end
