function M = stima3(vertices)
%STIMA3   Computes element stiffness matrix for simplex.
%   M = STIMA3(X) computes element stiffness matrix for simplex, 
%   i.e. for triangles in two dimensions (d=2) and tetraeder in 
%   three dimensions (d=3). The coordinates of the vertices are stored 
%   in X. X has dimension (d+1) x d. In two-dimension, the vertices 
%   are numbered anti-clockwise. M has dimension (d+1) x (d+1). In 
%   three-dimension, the vertices are numbered s.t. max(eig(M)) > 0. 
%
%   This routine should not be modified.

%    J. Alberty, C. Carstensen and S. A. Funken  02-11-99
%    File <stima3.m> in $(HOME)/acf/fem2d/ and
%                    in $(HOME)/acf/fem3d/ and
%                    in $(HOME)/acf/fem2d_heat/


%   Add psi and phi into stiffness matrix. Depends on psi and phi's
%   resolution. 

%   psi base(i)_x base(j)_x + phi base(i)_y base(j)_y
%   derivatives of base(i) base(j) are constants. thus only have to find
%   the average of function psi and phi on each elements.
%   
%   Integration on a trianle element for psi or phi. Use interpolation method is
%   preferred right now.


% first order
%   Int(phi)(a  ,b ,c) = 1/3(phi(a) + phi(b) + phi(c)) *Area. 
%
% third order
%   Int(phi)(a  ,b ,c) = 27/60 phi((a+b+c)/3) + 1/20(phi(a) + phi(b) +
%   phi(c)) + 2/15 (phi((a+b)/2) + phi((b+c)/2) + phi((c+a)/2) )

d = size(vertices,2);
G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
M = det([ones(1,d+1);vertices']) * G * G' / prod(1:d);
