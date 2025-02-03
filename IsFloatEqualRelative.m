% ----------------------------------------------------------------------------------
% Function for comparing (normalized) numbers with relative error maxRelDiff
%
% Input
% a,b - compared floating point numbers
% maxRelDiff - maximum permissible relative error
% Output:
% flag = 1/0 - numbers are equal/not equal
% --------------------------------------------------------------------------
function flag = IsFloatEqualRelative(a, b, maxRelDiff)
    flag = abs(a - b) <= max([abs(a),abs(b)])*maxRelDiff;
end
