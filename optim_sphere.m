% -------------------------------------------------------------------- 
% Determining the optimal configuration of "n" points on a sphere
% --------------------------------------------------------------------
clc
close all
clear all

format short g

% Let's remove extra degrees of freedom of the system of points placed on the sphere 
% by fixing the sphere in some coordinate system.
% Let's fix the sphere with points in some coordinate system.
% The center of the sphere always coincides with the origin of the coordinate system. 
% Therefore, we assume that the sphere, like a rigid body, has 1 point (the center of the sphere)
% already fixed.
% The 2nd fixed point of the sphere is the 1st point placed on the sphere with angles: fi_(1)=0; teta_(1)=0;
% Now the sphere can rotate only around the OZ axis connecting the center and the 1st placed point.
% For final fixation of the sphere, we set the azimuth angle of the 2nd point to be constant: fi_(2)=0;

% Radius of a sphere
r=1;

% Number of points on a sphere
n = 5;

% STD of distance measurement
sigma = 0.01;
% STD of the sphere shape
sigma_m = 0.01;

% Number of angles involved in optimization
len_angels = 2*n - 3;
% The initial values of the angles are assigned randomly
angels = 2*pi*rand(1, len_angels);

% options=optimset('Display','iter','MaxFunEvals',10000,'MaxIter',10000,'TolX',1.0e-12,'TolFun',1.0e-12);
TolX = eps(r);
TolFun = eps(r);
options=optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',TolX,'TolFun',TolFun);
[angels,fval,exitflag] =  fminsearch(@(angels) optim_sphere_fun(angels, r, sigma, sigma_m),angels,options);
% angels
fval
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

% draw a sphere and a figure with an optimal configuration of points
paint(x,y,z,r);

% Matrix of distances between points
s = zeros(n,n);
for i=1:n-1
    for j=i+1:n
      s(i,j)=sqrt(((x(i) - x(j)).^2 + (y(i) - y(j)).^2 + (z(i) - z(j)).^2));
      s(j,i)=s(i,j);
    end 
end
disp('Matrix of distances between points'),disp(s);

% Matrix of half squares of distances between points
C = 0.5*s.^2;
disp('The sum of the squares of the distances from each point to all the others'), disp(sum(C));

% We perform pseudo-inversion of matrix C, leaving the 4 largest eigenvalues
[U,L]=eig(C);
di=diag(L);
[~, Ind]=sort(abs(di),'descend');
Ind1_4 = Ind(1:4);
b=ones(n,1);
bb=U'*b;
z=bb(Ind1_4)./di(Ind1_4);
V=U(:, Ind1_4);
x=V*z;
disp('Vector of solution of a system of linear equations'), disp(x');

% f=2*(x.^2)'*C*(x.^2);
% f
% b'*C*b

% Checking that the optimal tetrahedron has pairwise equal opposite edges
if n == 4
    c =sqrt(2*C);
    fprintf('Optimal tetrahedron has pairwise equal opposite edges\n');
    fprintf('c(1,2) - c(3,4)=%f\n',c(1,2) - c(3,4));
    fprintf('c(2,3) - c(1,4)=%f\n',c(2,3) - c(1,4));
    fprintf('c(1,3) - c(2,4)=%f\n',c(1,3) - c(2,4));
end
