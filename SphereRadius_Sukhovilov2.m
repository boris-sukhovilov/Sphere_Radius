function [roots, sigma_R, sigma_upper_bound, status] = SphereRadius_Sukhovilov2(R0, S, sigma, sigma_m, rmin0, rmax0, d)

    % Construction of the measured matrix of half-squared distances
    C = 0.5 * S.^2;
    n = size(C,1); 
    
    roots = find_roots(C, rmin0, rmax0, d);
    
    n_roots = length(roots);

    if n_roots == 0
       [root, ~, exitflag,~] = fzero(@(r) r_fun(r, n, C), R0);
%        [root, ~, exitflag,~] = fzero(@(r) r_fun(r,N,C), (rmin0 + rmax0)/2);
       if exitflag == 1
           roots = root;
           n_roots = 1;
       end
    end

%     n_roots = length(roots);

    if n_roots == 0
        status = 0;
        roots = NaN; sigma_R = NaN; sigma_upper_bound = NaN;
        return;
    else
        status = 1;
    end

   sigma_R = zeros(1, n_roots);
   sigma_upper_bound = zeros(1, n_roots);
   F = zeros(1, n_roots);
   for i = 1 : n_roots
       [sigma_R(i), sigma_upper_bound(i), F(i)] = rmse(roots(i), C, sigma, sigma_m);
   end
end

function f = r_fun(r, n, C)
    rr=(r.^2)*ones(n, n);
    L=eig(rr-C);
    [~,j]=sort(abs(L), 'descend'); 
    f=L(j(1))+L(j(2))+L(j(3))-n*r^2;
end

function [sigma_R, sigma_max, F] = rmse(r, C, sigma, sigma_m)
    % Calculating the standard deviation of the radius estimate

    n = size(C,1);

    [U,L]=eig((r .^2)*ones(n,n)-C);
    [~,ind]=sort(abs(diag(L)), 'descend');
    % let's take 3 vectors corresponding to the three largest eigenvalues
    U=U(:,[ind(1) ind(2) ind(3)]);
    
    summa=0;
    for i=1:3
      summa=summa+(sum(U(:,i))).^2;
    end 
    % summa

%     kov=zeros(3,3);
%     for i=1:3
%      for j=1:3
%       for k=1:n
%        for l=1:n
%         kov(i,j)=kov(i,j)+C(l,k)*(U(l,i)*U(k,i))*(U(l,j)*U(k,j));
%        end 
%       end
%      end
%     end

    kov=zeros(3,3);
    for i=1:3
     for j=1:3
        kov(i,j) =  sum(sum(C.*(U(:,i)*U(:,i)').*(U(:,j)*U(:,j)')));
     end
    end 
    
%     kov - kov1
    
    F=n-summa;
    % Definition of the standard deviation of the radius estimate - sigmar
    sigma_R=sqrt((sigma_m^2)./F+(2*sigma^2*(sum(sum(kov))))/(4*r^2*F^2));

    % Upper bound of the standard deviation of the radius estimate
    sigma_max=sqrt(sigma_m^2/F+18*(sigma/F)^2);

end

function roots = find_roots(C, rmin0, rmax0, d)

    N = size(C,1);

    % Calculating the number of isolation intervals
    if d == 0
        n = 0;
    else
        n = ceil((rmax0 - rmin0) / d);
    end
    
    % Initialization of a vector for storing roots
    roots = [];
    options = optimset('Display', 'off');

    for i = 1:n
        % Calculating the boundaries of the isolation interval
        x1 = rmin0 + (i-1) * d;
        x2 = rmin0 + i * d;

        % Checking the signs of a function at the boundaries of an interval
        if r_fun(x1,N,C) * r_fun(x2,N,C) <= 0
            % Call fzero to find the root in the interval [x1, x2]
            root = fzero(@(r) r_fun(r,N,C), [x1 x2], options);
            roots = [roots, root];
        end
    end

%     if isempty(roots)
%         fprintf('Roots not found in given interval\n');
%     else
%         fprintf('Roots:\n');
%         disp(roots);
%     end
end
