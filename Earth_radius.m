%  ----------------------------------------------------------- 
%     Determination Earth radius by flight distances
%  -----------------------------------------------------------
%
% https://www.travelmath.com/ (22.01.25)	
%               Moscow    Tokyo    Delhi     Montevideo  Anchorage  Sydney
% Moscow        0         7501     4341      13348       7017       14485
% Tokyo                   0        5850      18576       5572       7792
% Delhi                            0         15598       9175       10419
% Montevideo                                 0           13522      11882
% Anchorage                                              0          11802
% Sydney                                                            0

function Earth_radius()
    clc

    s = ...
       [0    7501  4341   13348   7017  14485
        7501  0    5850   18576   5572  7792
        4341  5850 0      15598   9175  10419
        13348 18576 15598 0       13522 11882
        7017  5572  9175  13522   0     11802
        14485 7792  10419 11882   11802 0];

    b = ones(size(s,1),1);

    % Initial approximation Earth radius
    r0=6000;
    
    % 1) Determining Earth radius by 4 points
    s4 = s(1:4,1:4);
    b4 = b(1:4,:);
    [r,~]=fzero(@(r) (1-b4'*((2.*(sin(s4./(2*r))).^2)\b4)), r0);
    fprintf('Earth radius by 4 points:%g\n', r);

    % 2) Determining Earth radius by 5 points
    s5 = s(1:5,1:5);
    b5 = b(1:5,:);
    [r,~]=fzero(@(r) (1-b5'*pseudo_inv(2.*(sin(s5./(2*r))).^2, 4)*b5), r0);
    fprintf('Earth radius by 5 points:%g\n', r);

    % 3) Determining Earth radius by 6 points
    [r, ~]=fzero(@(r) (1-b'*pseudo_inv(2.*(sin(s./(2*r))).^2, 4)*b), r0);
    fprintf('Earth radius by 6 points:%g\n', r);

    % Plot
    i=0;
    step=1;
    for r=r-10:step:r+10
        i=i+1;
        x(i)=r;
        y(i)=1-b'*pseudo_inv(2.*(sin(s./(2*r))).^2, 4)*b;
    end
    plot(x,y),grid
end

