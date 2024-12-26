function example_for_suh()
clc

sigma = 0.03/3;

sigma_m = 0.0;

format long

% ������ ��� ����������, ��� ����� ������� ����������� ����� ������������
% ����� ������ ��������� �������� ������� ������ �������
% example_for_suh

r=@(d)1/sqrt( sum(sum( inv(  [0 d(1:3); d(1) 0 d(4:5); d(2) d(4) 0 d(6);...
    d(3) d(5) d(6) 0].^2/2))));

% ��� ������ ��������
phi=0.0988; % ���� ���������� �� ������
% phi=0.2; % ���� ���������� �� ������
relerr=0.01*[1,1,0,0,0,1];     % ����� ������������� ����������� ���������

disp('������, ��� ����� ������ ��������� �������� ������� ������ �������')
disp(['���� ���������� �� ������  ',num2str(phi)])
disp('������������� ����������� ��� ���������')
disp(relerr)

% ����� ���������� ����� �� ����������� ����������
M1=[ sin(phi); 0;         cos(phi)]; % ��� ����� �������
M2=[-sin(phi); 0;         cos(phi)];  % ��� ����� ������
M3=[0;          sin(phi); -cos(phi)]; % �� ����� �� ���
M4=[0;         -sin(phi); -cos(phi)]; % �� ����� � ���
%  ���������� ����� ������� � ��� ������������
s=zeros(1,6);
s(1)=sqrt(     sum(   (M1-M2).^2   )         );%12
s(2)=sqrt(     sum(   (M1-M3).^2   )         );%13
s(3)=sqrt(     sum(   (M1-M4).^2   )         );%14
s(4)=sqrt(     sum(   (M2-M3).^2   )         );%23
s(5)=sqrt(     sum(   (M2-M4).^2   )         );%24
s(6)=sqrt(     sum(   (M3-M4).^2   )         );%34
disp('"������" ����������'); disp( num2str(s));

disp(['������ �� "������" ����������� (��������� �������):',num2str(r(s))])

% S = [0   s(1) s(2) s(3);
%     s(1) 0    s(4) s(5);
%     s(2) s(4) 0    s(6);
%     s(3) s(5) s(6) 0];

% [U,S,V] = svd(S)
% [V,D] = eig(S)


% disp(['������ �� "������" �����������   ',num2str(calculate_radius_Sukhovilov(S))])

s=s.*(1+relerr); % ��������� �������������

% s(2) = 0.5*(s(2)+s(5));
% s(5)=s(2);
% s(1) = 0.5*(s(1)+s(6));
% s(6)=s(1);

S = [0    s(1) s(2) s(3);
     s(1) 0    s(4) s(5);
     s(2) s(4) 0    s(6);
     s(3) s(5) s(6) 0];

tmp = (S(1,2) + S(3,4))/2;
% tmp = 0;
S(1,2) = tmp; 
S(3,4) = tmp; 
S(2,1) = tmp; 
S(4,3) = tmp; 

tmp = (S(1,4) + S(2,3))/2;
S(1,4) = tmp; 
S(2,3) = tmp; 
S(4,1) = tmp; 
S(3,2) = tmp; 

tmp = (S(1,3) + S(2,4))/2;
S(1,3) = tmp; 
S(2,4) = tmp; 
S(3,1) = tmp; 
S(4,2) = tmp; 

% [U,S,V] = svd(ss)
% [V,D] = eig(ss)

disp('���������� ����������'); disp( num2str(s));

% disp(['������ �� ���������� �����������   ',num2str(calculate_radius_Sukhovilov(ss))])

    [R1, x, sigma_R] = SphereRadius_Sukhovilov(S, sigma, sigma_m);
    fprintf('������ ��������� �����, ����������� 1-� �������: %g\n', R1);
    fprintf('RMSE of R: %g\n\n', sigma_R);
%     x'
%     1./sqrt(sum(x))

    rmin0 = R1 - 3*sigma_R; 
    rmax0 = R1 + 3*sigma_R;
    rmin0 = max([eps, rmin0]);
    rmin0, rmax0
    
    % ����������� ������ ��������� �������� �����
    d = 3*sqrt(sigma^2 + sigma_m^2)/1000;

    [R2, sigma_R, sigma_max, status] = SphereRadius_Sukhovilov2(R1, S, sigma, sigma_m, rmin0, rmax0, d);
    if status == 1
        for i = 1 : length(R2)
            fprintf('Metod 2 Radius: %g\n', R2(i));
            fprintf('RMSE of R: %g\n', sigma_R(i));
            fprintf('Upper bound for RMSE of R: %g\n', sigma_max(i));
            sm = sensitivity_matrix(R2(i), S);
            disp(sm);
        end
    else
        fprintf('Metod 2 Radius not found!\n');
    end

%     disp(['������ �� ���������� ����������� (��������� �������):',num2str(r(s))])

    % ������� ����� ���������
    [a, b, c, a1, b1, c1] = getTetrahedronEdges(S);
    
    % ���������� ������� ����� ������������ ����-�������: https://ru.wikipedia.org/wiki/��������
    R_Cayley_Menger = circumscribedSphereRadius_Cayley_Menger(a, b, c, a1 , b1, c1);
    fprintf('������ ��������� �����, ����������� ����� ������������ ����-�������: %g\n', R_Cayley_Menger);

    R_Grelle = circumscribedSphereRadius_Grelle(a, b, c, a1, b1, c1);
    fprintf('������ ��������� �����, ����������� �� ������� Grelle: %g\n', R_Grelle);

end
