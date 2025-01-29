function [sigma_lambda] = sigma_eigenvalues(S2, sigma, sigma_m, R0, fPrint)
   [V, L] = eig(S2);
   lambda = abs(diag(L));
   [~,ind] = sort(lambda, 'descend');
   
   sigma_lambda = zeros(1,4);
   for i=1:4
       k = ind(i);
       y = V(:,k).*V(:,k);
       lambda_i = lambda(k);
       disp_s = 2*sigma^2*y'*S2*y;
       disp_m = 4*sigma_m^2*(lambda_i/R0)^2*sum(y.*y);
       
       if fPrint == 1
           fprintf('Eigen value number:%d\tRMSE from error distance:%g\tRMSE from model error:%g\n', i, sqrt(disp_s), sqrt(disp_m));
       end
       
       sigma_lambda(i) = sqrt(disp_s+disp_m);
    end 

end