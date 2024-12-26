function search_roots()
    clc
    close all
    
%     R = 1;
%     numPoints = 4;
%     
%     h = R*0.001;   % ���������� �� ������ ����� �� �������������� ����������, ������������ �����
%     h
%     
%     theta = pi / 4; % ���� ����� �������������� �����������
% 
%     type_distrib = 1;
%     delta = 0.01;
%     delta_m = 0.0;
% %     delta_m = 0;
%     sigma_m = delta_m/3;
%     
%     thetaRange = [0 2*pi];
%     phiRange = [0 pi/2];
%     
%     rng('shuffle');
% 
%     generate_type = 0;
%     
%     if generate_type == 0
%         points = generateRandomPointsInSolidAngle(R, numPoints, thetaRange, phiRange, sigma_m);
%     else
%         points = generate_optim_tetrahedron_points(R, h, theta);
%     end
% 
%     % ��������� ������� ����������
%     [S, S0, sigma] = generateMatrixDistance(points, delta, type_distrib);
%     
%     [R_Sukhovilov1, x, sigma_R] = SphereRadius_Sukhovilov(S, sigma, sigma_m);
%     fprintf('������ ��������� �����, ����������� 1-� �������: %g\n', R_Sukhovilov1);
%     fprintf('RMSE of R: %g\n', sigma_R);
%     
%     rmin0 = R_Sukhovilov1 - 3*sigma_R; 
%     rmax0 = R_Sukhovilov1 + 3*sigma_R;
%     rmin0, rmax0
%     
% %     rmin0 = 0.9;
% %     rmax0 = 1.1;
%     
%     rmin0 = max([eps, rmin0]);

    r = 1;
    S = example_for_suh(r);   

    S2 = 0.5 * S.^2;
    n = size(S,1);
    
    rmin0 = 0.5; rmax0 = 1.2;
    step = 0.001;
    i = 1;
    for r = rmin0 : step : rmax0
        x(i) = r;
%         y(i) = r_det(r,n,S2);
        y(i) = r_fun(r,n,S2);
        i = i + 1;
    end
    figure
    plot(x,y); grid;
end

function f = r_det(r,n,C)
    rr=(r.^2)*ones(n,n);
    f = det(rr-C);
end

function f = r_fun(r,n,C)
    rr=(r.^2)*ones(n,n);
    L=eig(rr-C);   %L-����������� ��������
%     di=diag(L);
    [~,j]=sort(-abs(L));  % �-�� 'sort' ��������� � ������� �����������
                          % j-������ �������������� ������� ��������������� ���������
    f=L(j(1))+L(j(2))+L(j(3))-n*r^2;
end


function S = example_for_suh(r)

% ������ ��� ����������, ��� ����� ������� ����������� ����� ������������
% ����� ������ ��������� �������� ������� ������ �������
% example_for_suh

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

% x(1) = M1(1); x(2) = M2(1); x(3) = M3(1); x(4) = M4(1);
% y(1) = M1(2); y(2) = M2(2); y(3) = M3(2); y(4) = M4(2);
% z(1) = M1(3); z(2) = M2(3); z(3) = M3(3); z(4) = M4(3);

% r = 1;
% figure
% paint(x,y,z,r);

%  ���������� ����� ������� � ��� ������������
s=zeros(1,6);
s(1)=sqrt(     sum(   (M1-M2).^2   )         );%12
s(2)=sqrt(     sum(   (M1-M3).^2   )         );%13
s(3)=sqrt(     sum(   (M1-M4).^2   )         );%14
s(4)=sqrt(     sum(   (M2-M3).^2   )         );%23
s(5)=sqrt(     sum(   (M2-M4).^2   )         );%24
s(6)=sqrt(     sum(   (M3-M4).^2   )         );%34
disp('"������" ����������'); disp( num2str(s));

s=s.*(1+relerr); % ��������� � ������������

s(2) = 0.5*(s(2)+s(5));
s(5)=s(2);

disp('���������� c ������������� ���������'); disp( num2str(s));

S = [0    s(1) s(2) s(3);
     s(1) 0    s(4) s(5);
     s(2) s(4) 0    s(6);
     s(3) s(5) s(6) 0];

end
