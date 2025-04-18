% --------------------------------------------------------------------------------
% Calculating the radius of a sphere circumscribing a tetrahedron using 
% the Carnot method
% https://archive.org/details/jstor-2973351/page/n1/mode/2up
% --------------------------------------------------------------------------------
function R = SphereRadius_Carnot(a, b, c, a1, b1, c1)
    % a, b, c - lengths of the edges of the tetrahedron between vertices 1-2, 1-3, 2-3
    % a1, b1, c1 - lengths of the opposite edges of the tetrahedron between vertices 3-4, 2-4, 1-4  
    %          4
    %         /|\
    %        / | \
    %     c1/a1|  \b1
    %      /   3   \ 
    %     / b / \ c \
    %    1------------2
    %          a
         
    a = a^2;  a1 = a1^2;
    b = b^2;  b1 = b1^2;
    c = c^2;  c1 = c1^2;
    numerator = (a*a1)^2 + (b*b1)^2 + (c*c1)^2 - 2*(a*a1*b*b1+a*a1*c*c1 +b*b1*c*c1);
    denominator =4*(a1^2*a + a^2*a1 + b1^2*b + b^2*b1 + c1^2*c + c^2*c1 + ...
        a*b1*c1 + c*a1*b1 + b*a1*c1 + a*b*c - a*b*b1 - b*c*b1 - b*a1*b1 - a*a1*b1 - b*b1*c1 - ...
        c*b1*c1 - a*b*a1 - a*c*a1 - b*c*c1 - a*c*c1 - a*a1*c1 - c*a1*c1);
    
    R = sqrt(numerator / denominator);
end
