% -----------------------------
% Perform Monte Carlo tests
% -----------------------------
function [R_N_max, R_N_mean, R_N_sigma] = MonteKarloTest( ...
    NTEST, R, numPoints, sigma, sigma_m, R0, confidence_interval, ...
    phiRange, thetaRange, h, beta, generate_type, fPrint)

    for i = 1 : NTEST
        R_N = SphereRadiusFromDistance(R, numPoints, sigma, sigma_m, R0, confidence_interval, ...
              phiRange, thetaRange, h, beta, generate_type, fPrint);
        
        if i == 1
            R_N_max = zeros(size(R_N));
            R_N_mean = zeros(size(R_N));
            R_N_sigma = zeros(size(R_N));
        end
        
        for j = 1 : length(R_N_max)
            if isnan(R_N(j))
                R_N_max(j) = NaN;
                R_N_mean(j) = NaN;
                R_N_sigma(j) = NaN;
                continue;
            end
            d = abs(R-R_N(j));
            if d > R_N_max(j)
                R_N_max(j) = d;
            end
            R_N_mean(j) = R_N_mean(j) + d;
            R_N_sigma(j) = R_N_sigma(j) + d^2;
        end
    end
    
    for j = 1 : length(R_N_max)
        R_N_mean(j) = R_N_mean(j)/NTEST;
        if NTEST == 1
            R_N_sigma(j) = sqrt(R_N_sigma(j)/NTEST);
        else
            R_N_sigma(j) = sqrt(R_N_sigma(j)/(NTEST-1));
        end
    end
    
end
