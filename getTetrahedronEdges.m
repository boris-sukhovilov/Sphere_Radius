function [a, b, c, a1, b1, c1] = getTetrahedronEdges(S)
    % a, b, c - edges of the tetrahedron emanating from vertex 1
    a = S(1,2);
    b = S(1,3);
    c = S(1,4);
    % a1, b1, c1 - edges of the tetrahedron opposite to edges a, b, c
    a1 = S(3,4);
    b1 = S(2,4);
    c1 = S(2,3);
end
