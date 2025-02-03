% ----------------------------------------------------------------------------
% Calculation of optimal coordinates of n points on a sphere of radius r,
% providing the minimum root mean square deviation of the radius estimate
% 
% r - sphere radius
% n - number of points on a sphere
% points - matrix vectors of coordinates points on a sphere (size 3 x n)
% sigma - RMSE of distance measurement
% sigma_m - RMSE of the sphere shape
% ----------------------------------------------------------------------------
function [points, exitflag] = optimal_n_points_on_sphere_fix(r, n, sigma, sigma_m)
    % Let's remove extra degrees of freedom of the system of points placed on the sphere 
    % by fixing the sphere in some coordinate system.
    % Let's fix the sphere with points in some coordinate system.
    % The center of the sphere always coincides with the origin of the coordinate system. 
    % Therefore, we assume that the sphere, like a rigid body, has 1 point (the center of the sphere)
    % already fixed.
    % The 2nd fixed point of the sphere is the 1st point placed on the sphere with angles: fi_(1)=0; teta_(1)=0;
    % Now the sphere can rotate only around the OZ axis connecting the center and the 1st placed point.
    % For final fixation of the sphere, we set the azimuth angle of the 2nd point to be constant: fi_(2)=0;

    % Number of angles involved in optimization
    len_angels = 2*n - 3;
    % The initial values of the angles are assigned randomly
%     angels = 2*pi*rand(1, len_angels);

% % Fixing the sphere
% fi_(1)=0;
% teta_(1)=0;
% fi_(2)=0;
% teta_(2)=angels(1);

    alfa = pi/3;
     % polar angle 
    theta0 = acos(cos(alfa) + (cos(pi-alfa) - cos(alfa)) * rand);


    % teta_(2)=angels(1);
%     angels(1) = pi/2;
    angels(1) = theta0;
    % fi_(3)=0;
    angels(2) = 2*pi/3;
    % teta_(3)=0;
%     angels(3) = pi/2;
    angels(3) = theta0;
    % fi_(4)=0;
    angels(4) = 4*pi/3;
    % teta_(4)=angels(1);
%     angels(5) = pi/2;
    angels(5) = theta0;
    
    num_rand_generate = (len_angels - 5)/2;
    
    % Generate random angles in given ranges
    phiRange = [0 2*pi];
     % azimuth angle
    phi = phiRange(1) + (phiRange(2) - phiRange(1)) * rand(num_rand_generate, 1);
    thetaRange = [0 pi];
     % polar angle 
    theta = acos(cos(thetaRange(1)) + (cos(thetaRange(2)) - cos(thetaRange(1))) * rand(num_rand_generate, 1));
   
    j = 5;
    for i = 1 : num_rand_generate
        j = j + 1;
        % fi angle
        angels(j) = phi(i);
        % teta angle
        j = j +1;
        angels(j) = theta(i);
    end
    

    TolX = eps(r);
    TolFun = eps(r);
    options=optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',TolX,'TolFun',TolFun);
    [angels,fval,exitflag] =  fminsearch(@(angels) optim_sphere_fun(angels, r, sigma, sigma_m),angels,options);
    % angels
    % fval
    % exitflag

    % Fixing the sphere
    % fixing the 1st point
    fi_(1)=0;
    teta_(1)=0;
    % partially fixing the 2nd point
    fi_(2)=0;
    teta_(2)=angels(1);

    j=1;
    for i = 3 : n
        j=j+1;
        fi_(i)=angels(j);
        j=j+1;
        teta_(i)=angels(j);
    end

    % Vectors of coordinates of points on a sphere
    x = zeros(1,n);
    y = zeros(1,n);
    z = zeros(1,n);
    for i=1:n
        x(i)=r*sin(teta_(i))*cos(fi_(i));
        y(i)=r*sin(teta_(i))*sin(fi_(i));
        z(i)=r*cos(teta_(i));
    end
 
    points = [x; y; z];
end
