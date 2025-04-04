% -------------------------------------------------------------------------------------
% Calculating the radius of a sphere circumscribing a tetrahedron using 
% the Euler, Grelle method
% https://archive.org/details/jstor-2973351/page/n1/mode/2up
% In 1752 Euler gave in effect, the following expression for V
% We use this relation in form Killing and Hovestadt state 
% In 1821 Crelle was published formula for R
% -------------------------------------------------------------------------------------
function R = SphereRadius_Euler_Grelle(a, b, c, a1, b1, c1)
    % a, b, c - lengths of the edges of the tetrahedron between vertices 1-2, 1-3, 1-4
    % a1, b1, c1 - lengths of the opposite edges of the tetrahedron between vertices 3-4, 2-4, 2-3
    %          4
    %         /|\
    %        / | \
    %     c /a1|  \b1
    %      /   3   \ 
    %     / b / \c1 \
    %    1-----------2
    %          a    

    V = (1/12)*sqrt((a^2+b^2+c^2+a1^2+b1^2+c1^2)*(a^2*a1^2+b^2*b1^2+c^2*c1^2) - ...
        2*a^2*a1^2*(a^2+a1^2) - 2*b^2*b1^2*(b^2+b1^2) - 2*c^2*c1^2*(c^2+c1^2) - ...
        a1^2*b^2*c^2 - a^2*b1^2*c^2 - a^2*b^2*c1^2 - a1^2*b1^2*c1^2);
    
    p = (a*a1+b*b1+c*c1)/2;
    R = sqrt(p*(p-a*a1)*(p-b*b1)*(p-c*c1))/(6*V);
end
