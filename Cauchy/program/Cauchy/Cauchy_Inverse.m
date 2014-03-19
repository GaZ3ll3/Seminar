function [u] = Cauchy_Inverse(geom, hmax)
% Solve Inverse Cauchy Problem on any geometry mesh. Default mesh is unit
% square.

% Now the problem is
% \Delta u = f 
%        u = u_d      on Dirichlet BC
%    du/dn = g        on Neumann BC

% This is an inverse solver.

global A b coordinates FreeNodes dirichlet neumann Ddata dset dnum elements3 AREA

if nargin < 1
    % default setting is UnitSquare
    % This is for modifying, Problem on 1 x r rectanle.
    r = 1.0;
    geom = [2 0 1 0 0 1 0;
        2 1 1 0 r 1 0;
        2 1 0 r r 1 0;
        2 0 0 r 0 1 0]';

    hmax = .02;
end


[A,b,coordinates,elements3, dirichlet, neumann, FreeNodes,AREA] = Cauchy_Init(geom,hmax);

% Dirichlet condition to be matched for Neumann boundary also as exact
% solution

% sphere
% Ddata = 0.500*(coordinates(:,1).^2 + coordinates(:,2).^2); 

% plane
% Ddata = coordinates(:,1) + coordinates(:,2);

% hyperbolic sine product cosine surface
% Ddata = sinh(coordinates(:,1)).*cos(coordinates(:,2)); 



% TEST DATA
% Ddata =  coordinates(:,2) + coordinates(:,2).^2/4;

g1 = 1;
g2 = 1;
g =sqrt(g1^2 + g2^2);
V0 = .5;
Ddata =  coordinates(:,1) + (V0+ g1*coordinates(:,1))*g1.*(1 - cosh(g*coordinates(:,2)))./(g*(g*cosh(g*coordinates(:,2)) ...
    - g2*sinh(g*coordinates(:,2))));

%-------------------------------------------------------------------------

% Use 'MaxFunEvals' iff 'GradObj' is off.
% options = optimset('Diagnostics','off','DerivativeCheck','off','FinDiffType','central','LargeScale','off',...
%     'GradObj','on','Display','iter-detailed','TolFun',1e-24,'TolX',1e-24,'MaxIter',10000,'MaxFunEvals',1000000,'HessUpdate','bfgs');


options = optimset('GradObj','on','Hessian','on','Display','iter-detailed','TolFun',1e-16,'TolX',1e-16,'MaxIter',10000,'MaxFunEvals',1000000);
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


subplot(1,3,1);
    trisurf(elements3,coordinates(:,1),coordinates(:,2),0*u',u','edgecolor','none','facecolor','interp');
    view(2);daspect([1 1 1 ]);
    title('Numerical Solution');
    colorbar;
subplot(1,3,2);
    trisurf(elements3,coordinates(:,1),coordinates(:,2),0*Ddata',Ddata','edgecolor','none','facecolor','interp');
    view(2);daspect([1 1 1]);
    title('Exact Solution');
    colorbar;
subplot(1,3,3);
    trisurf(elements3,coordinates(:,1),coordinates(:,2),0*Ddata',Ddata' - u','edgecolor','none','facecolor','interp');
    view(2);daspect([1 1 1]);
    title('Abs Error');
    colorbar;

save('numerical.mat','u');save('exact.mat','Ddata');

ERR = .0;
TRU = .0;
level = 1.0;

for i = 1 : size(elements3,1)
    if (coordinates(elements3(i,1),2) <= level) && ((coordinates(elements3(i,2),2) <= level)) && (coordinates(elements3(i,3),2) <= level) 
        
        ERR = ERR + AREA(i)*((Ddata(elements3(i,1)) - u(elements3(i,1)))^2 + ...
            (Ddata(elements3(i,2)) - u(elements3(i,2)))^2 + ...
            (Ddata(elements3(i,3)) - u(elements3(i,3)))^2 + ... 
            abs(Ddata(elements3(i,1)) - u(elements3(i,1)))*abs(Ddata(elements3(i,2)) - u(elements3(i,2))) +...
            abs(Ddata(elements3(i,2)) - u(elements3(i,2)))*abs(Ddata(elements3(i,3)) - u(elements3(i,3))) +...
            abs(Ddata(elements3(i,3)) - u(elements3(i,3)))*abs(Ddata(elements3(i,1)) - u(elements3(i,1))))/6;
        TRU = TRU + AREA(i)*((Ddata(elements3(i,1)))^2 +...
                               (Ddata(elements3(i,2)))^2 + ...
                               (Ddata(elements3(i,3)))^2 + ...
                               abs(Ddata(elements3(i,1)))* abs(Ddata(elements3(i,2)))  +...
                               abs(Ddata(elements3(i,2)))* abs(Ddata(elements3(i,3))) + ...
                               abs(Ddata(elements3(i,3)))* abs(Ddata(elements3(i,1))))/6;
    end
end

disp('Relatvie L2 error of solution and exact solution on unitsquare is:');
disp(sqrt(ERR)/sqrt(TRU));
end


function [f, g, Hess] = Cauchy_Opt(u_d)


global A b coordinates FreeNodes dirichlet neumann Ddata dset dnum elements3

% Regularization. 1e-11 or 1e-12 at most right now. will give an up-to 25%
% related abs error.
eps = 1e-12;
gleps = 1e-11;

f = .0; f_reg = .0; 
g = zeros(size(u_d,1),1);
g_reg = sparse(size(A,1),size(dirichlet,1)); 
% points x edges matrix
reg = sparse(size(dirichlet,1),1); % col base
delta = sparse(size(dirichlet,1),1);
dreg = sparse(size(dirichlet,1),1);
tg_reg = sparse(size(dirichlet,1),size(A,1));
% edges x points matrix

BC = sparse(size(A,1),1);
BC(dset) = u_d;

[~,~,u] = Cauchy_Reduction(A,b,coordinates,FreeNodes,dirichlet,neumann, BC);

figure(1);
trisurf(elements3,coordinates(:,1),coordinates(:,2),0*u',u','edgecolor','none','facecolor','interp');
view(2);daspect([1 1 1]);
colorbar;title('Numerical Solution of the Problem');


% UNCOMMENT
for j = 1 : size(neumann,1)
  f = f + norm(coordinates(neumann(j,1),:)- ...
      coordinates(neumann(j,2),:)) * ((u(neumann(j,1)) - Ddata(neumann(j,1)))^2 + (u(neumann(j,2)) - Ddata(neumann(j,2)))^2)/2/2;
end


% Regularizations



for j = 1 : size(dirichlet,1)
    delta(j) = norm(coordinates(dirichlet(j,1),:) - coordinates(dirichlet(j,2),:));
    reg(j) = (u(dirichlet(j,1)) - u(dirichlet(j,2)))/delta(j);
    dreg(j) = reg(j)/delta(j);
    f_reg = f_reg + (reg(j))^2;
end

%UNCOMMENT
f = f + eps*f_reg;

f = f + gleps*u'*A*u; % more regularization




Derivative = sparse(size(A,1), size(A,1));

% matrix of n x d



Derivative(dset,dset) = eye(dnum);
Derivative(FreeNodes, dset) = -A(FreeNodes, FreeNodes)\A(FreeNodes,dset);

Hess = zeros(dnum,dnum);

% UNCOMMENT
for j = 1:size(neumann,1)
    
    g = g + norm(coordinates(neumann(j,1),:)- ...
      coordinates(neumann(j,2),:)) * ((u(neumann(j,1)) - Ddata(neumann(j,1)))*Derivative(neumann(j,1),dset)' +...
             (u(neumann(j,2)) - Ddata(neumann(j,2)))*Derivative(neumann(j,2),dset)'    )/2      ;
        
    Hess = Hess + norm(coordinates(neumann(j,1),:)- ...
      coordinates(neumann(j,2),:)) *(Derivative(neumann(j,1),dset)'*Derivative(neumann(j,1),dset) +...
      Derivative(neumann(j,2),dset)'*Derivative(neumann(j,2),dset))/2;
    
end 
 

for j = 1 : size(dirichlet,1)
    g_reg(dirichlet(j,1),j) = 1;
    g_reg(dirichlet(j,2),j) = -1;
    tg_reg(j,dirichlet(j,1)) = 1/delta(j);
    tg_reg(j,dirichlet(j,2)) = -1/delta(j);
end

g = g + 2*eps*g_reg(dset,:)*dreg;

g = g + 2*gleps*Derivative(:,dset)'*A*u;

Hess = Hess + 2*eps*(tg_reg(:,dset)'*tg_reg(:,dset));

Hess = Hess + 2*gleps*Derivative(:,dset)'*A*Derivative(:,dset);




end

