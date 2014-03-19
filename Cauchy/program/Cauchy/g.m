function Stress = g(x)
%G   Data on the Neumann boundary
% U = 0.5*(x^2 + y^2)
Stress = zeros(size(x,1),1);

% U = x + y
% Stress = -ones(size(x,1),1);

% U = sinh(x)*cos(y)
% Stress = sinh(x(1)).*sin(x(2));
