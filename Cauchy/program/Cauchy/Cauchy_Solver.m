function [u] = Cauchy_Solver(geom,hmax)
% Split into two parts. One is settings, the other is solver.

if nargin < 1
    % default setting is UnitSquare
    geom = [ 2 0 1 0 0 1 0;
        2 1 1 0 1 1 0;
        2 1 0 1 1 1 0;
        2 0 0 1 0 1 0]';
    
    hmax = 0.02;
end
[A,b,coordinates,elements3, dirichlet, neumann, FreeNodes] = Cauchy_Init(geom,hmax);


u_d = 0.500*(coordinates(:,1).^2 + coordinates(:,2).^2); % boundary condition on other sides

% elements for showing

tic;
[~,~,u] = Cauchy_Reduction(A,b,coordinates,FreeNodes,dirichlet,neumann,u_d);
toc;


% Derivative = sparse(size(A,1), size(dirichlet,1));
% 
% % matrix of n x d
% 
% 
% 
% Derivative(unique(dirichlet),unique(dirichlet)) = eye(size(unique(dirichlet),1));
% Derivative(FreeNodes, unique(dirichlet)) = -A(FreeNodes, FreeNodes)\A(FreeNodes,unique(dirichlet));
% 
% grad = Derivative(unique(neumann),unique(dirichlet));

plt_flag = 1;

if (plt_flag == 1)
trisurf(elements3,coordinates(:,1),coordinates(:,2),u','facecolor','interp')
view(10,40);
title('Solution of the Problem');
end

end
