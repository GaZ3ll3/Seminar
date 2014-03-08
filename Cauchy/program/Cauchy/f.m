function VolumeForce = f(x)
%F   Volume force in considered domain.
% U = 0.5*(x^2 + y^2)
VolumeForce = -2*ones(size(x,1),1);

