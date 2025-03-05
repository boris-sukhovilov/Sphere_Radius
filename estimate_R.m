function R_hat = estimate_R(S, sigma)
    n = size(S, 1);
    var_r = compute_var_r(S, sigma);
    r = sqrt(mean(0.5 * S.^2, 2));  
    w = 1 ./ var_r;  
    R_hat = sum(w .* r) / sum(w);  
end
