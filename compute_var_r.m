function var_r = compute_var_r(S, sigma)
    n = size(S, 1);
    S2 = 0.5 * S.^2;  
    X = mean(S2, 2);  
    d_sq_mean = mean(S2, 2);  
    var_X = (sigma^2 / n^2) * sum(S.^2, 2) + (sigma^4 / (2 * n));
    var_r = var_X ./ (4 * d_sq_mean);  
end
