function [A,b,u] = Cauchy_Reduction(A,b,coordinates,FreeNodes,dirichlet,neumann, u_d)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for j = 1 : size(neumann,1)
    vec = coordinates(neumann(j,1),:) - coordinates(neumann(j,2),:);
%     b(neumann(j,:))=b(neumann(j,:)) + ...
%         (Psi(sum(coordinates(neumann(j,:),:))/2)*vec(2)^2 + Phi(sum(coordinates(neumann(j,:),:))/2)*vec(1)^2)...
%         *g(sum(coordinates(neumann(j,:),:))/2)/2/norm(vec);


% for speed up
%     b(neumann(j,:))=b(neumann(j,:)) + ...
%     Phi(sum(coordinates(neumann(j,:),:))/2)*norm(vec)...
%     *g(sum(coordinates(neumann(j,:),:))/2)/2;

% for accuracy and speed up
    x1 = coordinates(neumann(j,1),:);
    x2 = coordinates(neumann(j,2),:);
    theta = sqrt(3/5);
    mid = (x1+x2)/2;
    left = mid + (x2 - x1)*theta;
    right = mid - (x2 - x1)*theta;
    b(neumann(j,:)) = b(neumann(j,:)) + ...
        (4*Phi(mid)*g(mid)/9 + ...
        5*Phi(left)*g(left)/18 + ...
        5*Phi(right)*g(right)/18)*norm(vec);


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

