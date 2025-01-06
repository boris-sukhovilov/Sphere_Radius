function [roots, sigma_R, sigma_max, status] = SphereRadius_Sukhovilov2(R1, S, sigma, sigma_m, rmin0, rmax0, d)

    % 2-� ����� ������ ������� c �������� ������,
    % � �������� ��������� ���������� (������ ����� ����������� ��������)
    % ���������� ������� S2
    C = 0.5 * S.^2;
    N = size(C,1); 
    
%     [rmin, rmax, status] = find_isolation_interval(C, rmin0, rmax0);
%     fprintf('Isolation interval are:[%g  %g]\n', rmin, rmax);
    
%     if status == 0
%         r = NaN; sigma_R = NaN; sigma_max = NaN;
%         return;
%     end
%     
%     options = optimset('Display', 'on');
%     r = fzero(@(r) r_fun(r,N,C), [rmin rmax], options);
    
    roots = find_roots(C, rmin0, rmax0, d);
    
    n_roots = length(roots);

    if n_roots == 0
%        [root, ~, exitflag,~] = fzero(@(r) r_fun(r,N,C), R1);
       [root, ~, exitflag,~] = fzero(@(r) r_fun(r,N,C), (rmin0 + rmax0)/2);
       if exitflag == 1
           roots = root;
           n_roots = 1;
       end
    end

%     n_roots = length(roots);

    if n_roots == 0
        status = 0;
        roots = NaN; sigma_R = NaN; sigma_max = NaN;
        return;
    else
        status = 1;
    end

   sigma_R = zeros(1, n_roots);
   sigma_max = zeros(1, n_roots);
   F = zeros(1, n_roots);
   for i = 1 : n_roots
       [sigma_R(i), sigma_max(i), F(i)] = rmse(roots(i), C, sigma, sigma_m);
   end
end

function f = r_fun(r,n,C)
    rr=(r.^2)*ones(n,n);
    L=eig(rr-C);   %L-����������� ��������
%     di=diag(L);
    [~,j]=sort(-abs(L));  % �-�� 'sort' ��������� � ������� �����������
                          % j-������ �������������� ������� ��������������� ���������
    f=L(j(1))+L(j(2))+L(j(3))-n*r^2;
end

% function f = r_fun(r2,n,C)
%     rr=(r2)*ones(n,n);
%     L=eig(rr-C);   %L-����������� ��������
% %     di=diag(L);
%     [~,j]=sort(-abs(L));  % �-�� 'sort' ��������� � ������� �����������
%                            % j-������ �������������� ������� ��������������� ���������
%     f=L(j(1))+L(j(2))+L(j(3))-n*r2;
% end


function [sigma_R, sigma_max, F] = rmse(r, C, sigma, sigma_m)

    n = size(C,1);

    % ���������� ��� ������ �������
    [U,L]=eig((r .^2)*ones(n,n)-C);   % U-����������� ������� � L-����������� ��������
    [~,ind]=sort(-abs(diag(L)));  % �-�� 'sort' ��������� � ������� �����������
                                % j-������ �������������� ������� ��������������� ���������
    % ������ 3 �������, ��������������� ���� ���������� ������������ ���������
    U=U(:,[ind(1) ind(2) ind(3)]);
    
%     [U,S,V]=svds( (r .^2)*ones(n,n)-C,3 );   % U-����������� ������� � L-����������� ��������
%     U,S,V
    
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
    
    % sumkov=2*sigma^2*(sum(sum(kov))); % ����� ��������� �������������� �������

    % ����������� ��� ������ �������  - sigmar
    F=n-summa;
    sigma_R=sqrt((sigma_m^2)./F+(2*sigma^2*(sum(sum(kov))))/(4*r^2*F^2));

    % ������� ������� ��� ������ �������
    sigma_max=sqrt(sigma_m^2/F+18*(sigma/F)^2);

end

% function [rmin, rmax, status] = find_isolation_interval(C, rmin0, rmax0)
% 
%     N = size(C,1);
%     
%     % ��������� ��������
%     rmin = rmin0;
%     rmax = rmax0;
% 
%     % ����������� ����� ������������
%     epsilon = eps(rmin);
% 
%     % ��������� ���������� �������������
%     n = 10;
% 
%     % ���� ��� ���������� �����
%     found = false;
% 
%     while ~found
%         % ���������� ��������� �� ������������
%         r = linspace(rmin, rmax, n);
% 
%         % ����� ��������� � ������� ������� ������� �� ��������
%         for i = 1:n-1
%             if r_fun(r(i),N,C) * r_fun(r(i+1),N,C) < 0
%                 rmin = r(i);
%                 rmax = r(i+1);
%                 found = true;
%                 status = 1;
%                 return;
%             end
%         end
% 
%         % ���������� ���������� �������������
%         n = n * 2;
% 
%         % �������� ����� ������������
%         if (rmax - rmin) / n < epsilon
%             fprintf('����� ������������ ����� ������ ���������� ����� ������� � ��������� ������. ������ ���������.\n');
%             status = 0;
%             return;
%         end
%     end
% 
%     if ~found
%         fprintf('������ �� ������ � �������� ���������.\n');
%         status = 0;
%     end
% end

function roots = find_roots(C, rmin0, rmax0, d)

    N = size(C,1);

    % ���������� ���������� ���������� ��������
    if d == 0
        n = 0;
    else
        n = ceil((rmax0 - rmin0) / d);
    end
    
    % ������������� ������� ��� �������� ������
    roots = [];
    options = optimset('Display', 'off');

    % ���� �� ���������� ��������
    for i = 1:n
        % ���������� ������ ��������� ��������
        x1 = rmin0 + (i-1) * d;
        x2 = rmin0 + i * d;

        % �������� ������ ������� �� �������� ���������
        if r_fun(x1,N,C) * r_fun(x2,N,C) <= 0
            % ����� fzero ��� ���������� ����� � ��������� [x1, x2]
            root = fzero(@(r) r_fun(r,N,C), [x1 x2], options);
            % ���������� ���������� �����
            roots = [roots, root];
        end
    end

    % ����� ����������
    if isempty(roots)
        fprintf('����� �� ������� � �������� ���������.\n');
    else
        fprintf('��������� �����:\n');
        disp(roots);
    end
end



