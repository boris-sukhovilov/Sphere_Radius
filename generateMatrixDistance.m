function [S, S0, sigma] = generateMatrixDistance(points, delta, type_distrib)

    numPoints = size(points, 2);
    
    S0 = zeros(numPoints,numPoints);
    S = zeros(numPoints,numPoints);
    
    % Добавление ошибок измерений
    if type_distrib == 0
    %    sigma=sqrt((delta.^2)./12);  % СКО для равномерного распределения
        sigma=delta/sqrt(3);  % СКО для равномерного распределения
    else
        sigma=delta/3;  % СКО для нормального распределения
    end
    
    for i = 1 : numPoints - 1
        for j = i + 1 : numPoints
            S0(i,j) = norm(points(:,i) - points(:,j));
            S0(j,i) = S0(i,j);
            if type_distrib == 0
                err=delta.*(2.*rand - 1);
            else
                err=sigma*randn;
            end
            S(i,j)= S0(i,j)+err;
            S(j,i)=S(i,j);
        end
    end
    
end
