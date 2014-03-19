function gamma = Psi(x)
% INPUT N X 2 ARRAY
% OUTPUT N X 1 ARRAY
% gamma = ones(size(x,1),1);

% gamma = 1 + 0.5*x(2); % 1+ 0.5t


g1 = 1;
g2 = 1.0;
g = sqrt(g1^2 + g2^2);
V0 = .5;

gamma = g*(V0 + g1*x(1))/(g*cosh(g*x(2)) - g2*sinh(g*x(2))); 


end

