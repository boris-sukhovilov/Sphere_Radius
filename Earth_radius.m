%  ----------------------------------------------------------- 
%   Determination Earth radius by flight distances
%  -----------------------------------------------------------
%
% https://www.travelmath.com/ (22.01.25)	
%         Seattle       Moscow  Tokyo   Delhi   Montevideo      Anchorage       Sydney
% Seattle 0             8397    7712    11332   11258           3639            12454
% Moscow                0       7501    4341    13348           7017            14485
% Tokyo                         0       5850    18576           5572            7792
% Delhi                                 0       15598           9175            10419
% Montevideo                                    0               13522           11882
% Anchorage                                                     0               11802
% Sydney                                                                        0

function Earth_radius()
    clc
    
    % STD of arc distance measurement (km)
    sigma = 10/3;
    % STD of linear distance measurement: sigma_s = sigma*cos(L(i,j)/(2*R))
    % When "2*R" is much greater than "L(i,j)" -> cos(L(i,j)/(2*R)) -> 1
    % and sigma_s -> sigma
    
    % STD of the sphere shape (km)
    sigma_m = 10/3;
    % Initial approximation Earth radius (km)
    R0 = 6000;

    % flight distances (arc distance measurement)
    s = ...
       [0      8397   7712    11332   11258   3639   12454
        0      0      7501    4341    13348   7017   14485
        0      0      0       5850    18576   5572   7792
        0      0      0       0       15598   9175   10419
        0      0      0       0       0       13522  11882
        0      0      0       0       0       0      11802
        0      0      0       0       0       0      0];
    
    s = s + s';
    
    n = size(s,1);
    b = ones(n,1);
    
    for i = 4 : n
        L = s(1:i,1:i);
        B = b(1:i,:);
        [R, ~]=fzero(@(r) (2-B'*pseudo_inv((sin(L./(2*r))).^2, 4)*B), R0);
        S = 2*R*sin(L./(2*R));
        sigma_R = STD_R(R, S, sigma, sigma_m);
        fprintf('Earth radius by %d points: %g\tSTD R: %g\n', i, R, sigma_R);
    end
    
end

% ----------------------------------------------------------------------
% sigma_R - STD of radius R
% ----------------------------------------------------------------------
function sigma_R = STD_R(R, S, sigma, sigma_m)
    S2 = 0.5*S.*S;
    S2_pinv = pseudo_inv(2.*R^2*(sin(S./(2*R))).^2, 4);
    b = ones(size(S2,1),1);
    x = S2_pinv*b;
    % the component of variance (1/R^2) caused by errors in distance measurement
    d=(2*sigma^2)*((x.^2)'*S2*(x.^2)); 
    % component of variance (1/R^2) due to model errors
    d=d+4*sigma_m^2*x'*x/R^2;
    % root mean square error (STD) R
    sigma_R=0.5*(R^3)*sqrt(d);
end
