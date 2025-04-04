function example_for_suh()
    % (c) Vasiliev U.S example for Sukhovilov,
    % that it is undesirable to locate points near the poles
    % small measurement errors cause a large radius error?

    clc
    
    % anonymous function
    r=@(d)1/sqrt( sum(sum(inv([0 d(1:3); d(1) 0 d(4:5); d(2) d(4) 0 d(6); d(3) d(5) d(6) 0].^2/2))));

    % radius sphere
    R = 1;
    % initial approximation radius sphere
    R0 = R + 0.5*R;

%     numPoints = 4;
    
    % For a normal distribution, the �4sigma confidence interval covers approximately 99.9937% of all possible values.
    % This means that the confidence level for the �4sigma interval is about 99.9937%.
    confidence_interval = 4;
    
    sigma = 0.01;
    sigma_m = 0.01;

    % initial calculation data
    phi=0.0988; % Angle of deviation from the pole
    relerr=0.01*[1,1,0,0,0,1];     % we set the relative measurement error

    disp('example that small measurement errors cause large radius error?')
    disp(['angle of deviation from the pole ',num2str(phi)])
    disp('relative errors for measurements')
    disp(relerr)

    % I set the coordinates of the points according to Sukhovilov's suggestion
    M1=[ sin(phi); 0;      cos(phi)];   % north pole to the right
    M2=[-sin(phi); 0;      cos(phi)];   % north pole to the left
    M3=[0;       sin(phi); -cos(phi)];  % south pole away from us
    M4=[0;      -sin(phi); -cos(phi)];  % south pole toward us
    % distances between points with machine error
    s=zeros(1,6);
    s(1)=sqrt(     sum(   (M1-M2).^2   )         );%12
    s(2)=sqrt(     sum(   (M1-M3).^2   )         );%13
    s(3)=sqrt(     sum(   (M1-M4).^2   )         );%14
    s(4)=sqrt(     sum(   (M2-M3).^2   )         );%23
    s(5)=sqrt(     sum(   (M2-M4).^2   )         );%24
    s(6)=sqrt(     sum(   (M3-M4).^2   )         );%34
    
    disp('"exact" distances'); disp( num2str(s));

    disp(['radius by "exact" distances (anonymous function):',num2str(r(s))])
    
    s=s.*(1+relerr); % measurement with error

    % s(2) = 0.5*(s(2)+s(5));
    % s(5)=s(2);
    % s(1) = 0.5*(s(1)+s(6));
    % s(6)=s(1);

    S = [0    s(1) s(2) s(3);
         s(1) 0    s(4) s(5);
         s(2) s(4) 0    s(6);
         s(3) s(5) s(6) 0];

    disp('measured distances'); disp( num2str(s));
    disp(['radius by "measured" distances (anonymous function):',num2str(r(s))])

    fPrint = 1;
    generate_type = 1;
    calc_Radius(S, sigma, sigma_m, R0, confidence_interval, fPrint, generate_type);
    
%     fprintf('Sphere radius for optimal placement of %d points= %g\n', numPoints, radius_optimal_n_points_on_sphere(S));
    
%     a = S(1,2); b = S(1,3); c = S(1,4);
%     R = (2^(1/2)*(a^2 + b^2 + c^2)^(1/2))/4;
%     a,b,c
%     R
end
