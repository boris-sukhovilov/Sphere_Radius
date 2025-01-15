function [sigma_lambda] = sigma_eigenvalues(S, sigma, sigma_m, R0)
   S2 = 0.5 * S.^2;
   [V, L] = eig(S2);
   [~,ind]=sort(abs(diag(L)), 'descend');
   
   sigma_lambda = zeros(1,4);
   for i=1:4
       k = ind(i);
       disp_s = 2*sigma^2*sum(sum(S2.*(V(:,k)*V(:,k)').*(V(:,k)*V(:,k)')));
       disp_m = 4*R0^2*sigma_m^2*sum(V(:,i))^2;
%        fprintf('Eigen value number:%d\tRMSE from error distance:%g\tRMSE from model error:%g\n', i, sqrt(disp_s), sqrt(disp_m));
       sigma_lambda(i) = sqrt(disp_s+disp_m);
    end 

end