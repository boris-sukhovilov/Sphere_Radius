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
function points = optimal_n_points_on_sphere(r, n, sigma, sigma_m)
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
    angels = 2*pi*rand(1, len_angels);

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
