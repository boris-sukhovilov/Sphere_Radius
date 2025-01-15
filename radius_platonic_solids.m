function R = radius_platonic_solids(S)
    numPoints = size(S,1);
    S2 = 0.5*S.*S;
    R = sqrt(sum(sum(S2)))/numPoints;
end
