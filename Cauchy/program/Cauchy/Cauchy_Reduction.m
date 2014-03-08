function [A,b,u] = Cauchy_Reduction(A,b,coordinates,FreeNodes,dirichlet,neumann, u_d)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for j = 1 : size(neumann,1)
  b(neumann(j,:))=b(neumann(j,:)) + norm(coordinates(neumann(j,1),:)- ...
      coordinates(neumann(j,2),:)) * g(sum(coordinates(neumann(j,:),:))/2)/2;
end

% Dirichlet conditions 
u = sparse(size(coordinates,1),1);
% u(unique(dirichlet)) = u_d(coordinates(unique(dirichlet),:));
u(unique(dirichlet)) = u_d(unique(dirichlet));
b = b - A * u;

% Computation of the solution
u(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);

% graphic representation

% elements4 = [];
% show(elements3,elements4,coordinates,full(u));


end

