function [u] = Cauchy_Inverse(geom, hmax)
% Solve Inverse Cauchy Problem on any geometry mesh. Default mesh is unit
% square.

% Now the problem is
% \Delta u = f 
%        u = u_d      on Dirichlet BC
%    du/dn = g        on Neumann BC

% This is an inverse solver.

global A b coordinates FreeNodes dirichlet neumann Ddata dset dnum elements3

% Enable figure flags to visualize solutions.
plt_flag = 1; % plot numerical solution
plt_exact = 1; % plot exact solution
cmp = 1; % plot difference


if nargin < 1
    % default setting is UnitSquare
    geom = [ 2 0 1 0 0 1 0;
        2 1 1 0 1 1 0;
        2 1 0 1 1 1 0;
        2 0 0 1 0 1 0]';
    
    hmax = 0.03;
end


[A,b,coordinates,elements3, dirichlet, neumann, FreeNodes] = Cauchy_Init(geom,hmax);

% Dirichlet condition to be matched for Neumann boundary also as exact
% solution
Ddata = 0.500*(coordinates(:,1).^2 + coordinates(:,2).^2); 

% Use 'MaxFunEvals' iff 'GradObj' is off.
options = optimset('Diagnostics','off','DerivativeCheck','off','FinDiffType','central','LargeScale','off',...
    'GradObj','on','Display','iter-detailed','TolFun',1e-16,'TolX',1e-16,'MaxIter',10000,'MaxFunEvals',1000000,'HessUpdate','bfgs');
% initial guess

dset = unique(dirichlet);
dnum = size(dset,1);

% Initial guess, since inverse problem has been stablized. We can use rand
% function to generate data for Dirichlet BC.

u_d = rand(dnum,1);

% UNCOMMENT BELOW to use exact data for inversion process.
% u_d = Ddata(dset);

[u_d,~] = fminunc(@Cauchy_Opt,u_d,options);

UD = zeros(size(A,1),1);
UD(dset) = u_d;

[~,~,u] = Cauchy_Reduction(A,b,coordinates,FreeNodes,dirichlet,neumann,UD);




if (plt_flag == 1)
    figure(1);
    trisurf(elements3,coordinates(:,1),coordinates(:,2),u','facecolor','interp')
    view(10,40);
    title('Numerical Solution of the Problem');
end

if (plt_exact == 1)
    figure(2);
    trisurf(elements3,coordinates(:,1),coordinates(:,2),Ddata','facecolor','interp')
    view(10,40);
    title('Exact Solution of the Problem');
end

if (cmp == 1)
    figure(3);
    trisurf(elements3,coordinates(:,1),coordinates(:,2),Ddata' - u','facecolor','interp')
    view(10,40);
    title('Absolute Difference of solutions');
end
end


function [f, g] = Cauchy_Opt(u_d)


global A b coordinates FreeNodes dirichlet neumann Ddata dset dnum elements3

% Regularization.
eps = 1e-11;

f = .0; f_reg = .0; 
g = zeros(size(u_d,1),1);
g_reg = sparse(size(A,1),size(dirichlet,1)); % comprehension needed
reg = sparse(size(dirichlet,1),1); % col base
dreg = sparse(size(dirichlet,1),1);
% edges x points matrix

BC = sparse(size(A,1),1);
BC(dset) = u_d;

[~,~,u] = Cauchy_Reduction(A,b,coordinates,FreeNodes,dirichlet,neumann, BC);

figure(1);
trisurf(elements3,coordinates(:,1),coordinates(:,2),u','facecolor','interp')


for j = 1 : size(neumann,1)
  f = f + norm(coordinates(neumann(j,1),:)- ...
      coordinates(neumann(j,2),:)) * ((u(neumann(j,1)) - Ddata(neumann(j,1)))^2 + (u(neumann(j,2)) - Ddata(neumann(j,2)))^2)/2/2;
end


% Regularizations



for j = 1 : size(dirichlet,1)
    reg(j) = (u(dirichlet(j,1)) - u(dirichlet(j,2)))/norm(coordinates(dirichlet(j,1),:) - coordinates(dirichlet(j,2),:));
    dreg(j) = reg(j)/norm(coordinates(dirichlet(j,1),:) - coordinates(dirichlet(j,2),:));
    f_reg = f_reg + (reg(j))^2;
end

f = f + eps*f_reg;


Derivative = sparse(size(A,1), size(A,1));

% matrix of n x d



Derivative(dset,dset) = eye(dnum);
Derivative(FreeNodes, dset) = -A(FreeNodes, FreeNodes)\A(FreeNodes,dset);



% grad = Derivative(unique(neumann),unique(dirichlet));

for j = 1:size(neumann,1)
    
    g = g + norm(coordinates(neumann(j,1),:)- ...
      coordinates(neumann(j,2),:)) * ((u(neumann(j,1)) - Ddata(neumann(j,1)))*Derivative(neumann(j,1),dset)' +...
             (u(neumann(j,2)) - Ddata(neumann(j,2)))*Derivative(neumann(j,2),dset)'    )/2      ;
    
end


for j = 1 : size(dirichlet,1)
    g_reg(dirichlet(j,1),j) = 1;
    g_reg(dirichlet(j,2),j) = -1;
end

g = g + 2*eps*g_reg(dset,:)*dreg;




end

