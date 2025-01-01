function search_roots()
    clc
    close all
    
%     R = 1;
%     numPoints = 4;
%     
%     h = R*0.001;   % Расстояние от центра сферы до горизонтальных плоскостей, расположения точек
%     h
%     
%     theta = pi / 4; % Угол между диаметральными плоскостями
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
%     % Генерация матрицы расстояний
%     [S, S0, sigma] = generateMatrixDistance(points, delta, type_distrib);
%     
%     [R_Sukhovilov1, x, sigma_R] = SphereRadius_Sukhovilov(S, sigma, sigma_m);
%     fprintf('Радиус описанной сферы, вычисленный 1-м методом: %g\n', R_Sukhovilov1);
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
    L=eig(rr-C);   %L-собственные значения
%     di=diag(L);
    [~,j]=sort(-abs(L));  % ф-ия 'sort' сортирует в порядке возрастания
                          % j-вектор первоначальных номеров отсортированных элементов
    f=L(j(1))+L(j(2))+L(j(3))-n*r^2;
end


function S = example_for_suh(r)

% пример для Суховилова, что около полюсов располагать точки нежелательно
% малые ошибки измерений вызывают большую ошибку радиуса
% example_for_suh

% исх данные расчётов
phi=0.0988; % Угол отклонения от полюса
% phi=0.2; % Угол отклонения от полюса
relerr=0.01*[1,1,0,0,0,1];     % задаём относительную погрешность измерений

disp('пример, что малые ошибки измерений вызывают большую ошибку радиуса')
disp(['угол отклонения от полюса  ',num2str(phi)])
disp('относительные погрешности для измерений')
disp(relerr)

% задаю координаты точек по предложению Суховилова
M1=[ sin(phi); 0;         cos(phi)]; % сев полюс направо
M2=[-sin(phi); 0;         cos(phi)];  % сев полюс налево
M3=[0;          sin(phi); -cos(phi)]; % юж полюс от нас
M4=[0;         -sin(phi); -cos(phi)]; % юж полюс к нам

% x(1) = M1(1); x(2) = M2(1); x(3) = M3(1); x(4) = M4(1);
% y(1) = M1(2); y(2) = M2(2); y(3) = M3(2); y(4) = M4(2);
% z(1) = M1(3); z(2) = M2(3); z(3) = M3(3); z(4) = M4(3);

% r = 1;
% figure
% paint(x,y,z,r);

%  расстояния между точками с маш погрешностью
s=zeros(1,6);
s(1)=sqrt(     sum(   (M1-M2).^2   )         );%12
s(2)=sqrt(     sum(   (M1-M3).^2   )         );%13
s(3)=sqrt(     sum(   (M1-M4).^2   )         );%14
s(4)=sqrt(     sum(   (M2-M3).^2   )         );%23
s(5)=sqrt(     sum(   (M2-M4).^2   )         );%24
s(6)=sqrt(     sum(   (M3-M4).^2   )         );%34
disp('"точные" расстояния'); disp( num2str(s));

s=s.*(1+relerr); % измерения с погрешностью

s(2) = 0.5*(s(2)+s(5));
s(5)=s(2);

disp('расстояния c погрешностями измерения'); disp( num2str(s));

S = [0    s(1) s(2) s(3);
     s(1) 0    s(4) s(5);
     s(2) s(4) 0    s(6);
     s(3) s(5) s(6) 0];

end
