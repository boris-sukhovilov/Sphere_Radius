% // David Eberly, Geometric Tools, Redmond WA 98052
% // Copyright (c) 1998-2024
% // Distributed under the Boost Software License, Version 1.0.
% // https://www.boost.org/LICENSE_1_0.txt
% // https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
% // Version: 6.0.2023.08.08
% https://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf (5.1)
% https://www.geometrictools.com/GTE/Mathematics/ApprSphere3.h
% 
% Translate from C++ to Matlab by Copilot
function [center, radius, iterations] = FitUsingLengths(points, maxIterations, initialCenterIsAverage, epsilon)
    if nargin < 4
        epsilon = 0;
    end

    % Number of points
    numPoints = size(points, 1);

    % Calculating the average value of points
    average = mean(points, 1);

    % Initial guess for the center
    if initialCenterIsAverage
        center = average;
    else
        center = points(1, :);
    end

    epsilonSqr = epsilon^2;
    iterations = 0;

    for iteration = 1:maxIterations
        %Updating iterations
        current = center;

        % Calculation of average L, dL/da, dL/db, dL/dc
        lenAverage = 0;
        derLenAverage = zeros(1, 3);
        for i = 1:numPoints
            diff = points(i, :) - center;
            length = norm(diff);
            if length > 0
                lenAverage = lenAverage + length;
                invLength = 1 / length;
                derLenAverage = derLenAverage - invLength * diff;
            end
        end
        lenAverage = lenAverage / numPoints;
        derLenAverage = derLenAverage / numPoints;

        center = average + lenAverage * derLenAverage;
        radius = lenAverage;

        diff = center - current;
        diffSqrLen = dot(diff, diff);
        if diffSqrLen <= epsilonSqr
            break;
        end

        iterations = iterations + 1;
    end
end
