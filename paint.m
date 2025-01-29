% ****************************************************************
% Display the surface of a sphere with points on it
% Input: x,y,z - coordinates of points
% Output: none
% ****************************************************************
function paint(x,y,z,r)
% Display the surface of a sphere with radius r
[X,Y,Z] = sphere;
surf(r.*X,r.*Y,r.*Z)
axis equal;

[~, col]=size(x);

% Display points
hold on;
plot3(x,y,z,'*r','markersize',20);
delta = (1+.05);
for i = 1 : col
 text(delta*x(i), delta*y(i), delta*z(i), num2str(i),'Color','black','FontSize',16);
end

hold off;
