% *********************************************************************
% ���������� ���������� ������� � ��������� �� ������ �����������
% ��� ���������,�������������� ����� �������� �����
% *********************************************************************
clc

bb=[1 1 1 1]';
bb

syms a b c d e f real

% syms a12 a13 a14 a23 a24 a34 real
% s1=0.5*[0      a12^2   a13^2   a14^2
%        a12^2  0       a23^2   a24^2
%        a13^2  a23^2   0       a34^2
%        a14^2  a24^2   a34^2   0  ];


s=0.5*[0     a^2    b^2   c^2
       a^2   0      c^2   b^2
       b^2   c^2    0     a^2
       c^2   b^2    a^2   0  ];
   
eigen_values = eig(s);
simplify(eigen_values)

det(s)
simplify(det(s))

inv_s=inv(s);
inv_s

R=simplify(1/sqrt(bb'*inv_s*bb));
R

R2 = simplify(1/(bb'*inv_s*bb));
R2

x = simplify(inv_s*bb);
x

[V,D] = eig(s)

% eig(inv_s)

r = simplify(sqrt(1/sum(x)));
r
r*r

x2 = x.*x;
D = simplify(2*x2'*s*x2);
D

% D2 = 0;
% for i = 1 : 3
%     for j = i+1 : 4
%         D2 = D2 + 4*s(i,j)*x2(i)*x2(j);
%     end
% end
% D2 = simplify(D2);
% D2

D0 = 1./(8*(R^6));
D0

% ����������� (����������) ������, ����� ����� ����� � ���������
% �������
a = 1;
b = 1;
c = sqrt(2);

% �������������, ��������� ������ �������� 1
a = 1;
b = sqrt(4-1);
b
c = 2;

C=0.5*[0     a^2    b^2   c^2
       a^2   0      c^2   b^2
       b^2   c^2    0     a^2
       c^2   b^2    a^2   0  ];
   
inv_c = inv(C)
[U,S,V] = svd(C)

% R=1/sqrt(bb'*inv_c*bb);
R=1/sqrt(bb'*(C\bb));
R

% ����� ������� ����������
s=0.5*[0     a^2    b^2   c^2
       a^2   0      d^2   e^2
       b^2   d^2    0     f^2
       c^2   e^2    f^2   0  ];

inv_s=inv(s);
inv_s

bb=[1 1 1 1]';
R=simplify(1/sqrt(bb'*inv_s*bb));
R
