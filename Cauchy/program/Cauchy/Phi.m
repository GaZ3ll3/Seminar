function xi = Phi(x)
% INPUT N X 2 ARRAY
% OUTPUT N X 1 ARRAY
% xi = ones(size(x,1),1);


% xi = 1/(1 + 0.5*x(2));% 1 + 0.5t

g1 = 1;
g2 = 1.0;
g = sqrt(g1^2 + g2^2);
V0 = .5;

xi= (g*cosh(g*x(2)) - g2*sinh(g*x(2)))/(g*(V0 + g1*x(1))); 



end

