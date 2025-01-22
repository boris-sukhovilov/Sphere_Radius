% ------------------------------------------------------------------
% Minimized functional
%
% Input
% angels - vector of azimuthal and polar angles of points on the sphere
% r - radius of the sphere
% Output
% f - variance to be minimized
% ------------------------------------------------------------------
function [f] = optim_sphere_fun(angels, r)

% Number of points
n=(length(angels)+3)/2;

% Fixing the sphere
fi_(1)=0;
teta_(1)=0;
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

% Matrix of half squares of distances between points on a sphere
C = zeros(n,n);
for i=1:n-1
    for j=i+1:n
      C(i,j)=0.5*((x(i) - x(j)).^2 + (y(i) - y(j)).^2 + (z(i) - z(j)).^2);
      C(j,i)=C(i,j);
    end 
end 

% We perform a pseudo-inversion of the matrix C, leaving the 4 largest eigenvalues
[U,L]=eig(C);
di=diag(L);
[~, Ind]=sort(abs(di),'descend');
Ind1_4 = Ind(1:4);
b=ones(n,1);
bb=U'*b;
z=bb(Ind1_4)./di(Ind1_4);
V=U(:, Ind1_4);
Y=V*z;

% Functionality without taking into account deviations from the shape of a sphere
% Only measurement errors are taken into account
f=2*(Y.^2)'*C*(Y.^2);

% Functional taking into account the deviation from the shape of the sphere.
% The first term is determined by errors in measuring pairwise distances
% The second term is determined by errors in the shape of the surface of the sphere
% sigma = 1, sigma_m = 1
% f=2*(Y.^2)'*C*(Y.^2)+4/r^2*((Y.^2)'*b)^2;
    