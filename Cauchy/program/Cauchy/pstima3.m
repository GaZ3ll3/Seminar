function M = pstima3(vertices,cpsi,cphi)
% A COPY OF STIMA3, BUT WITH SUPPORT WITH PSI, PHI
d = size(vertices,2); % 2 dimensional
G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
Coef = diag([cpsi,cphi]);
M = det([ones(1,d+1);vertices']) * G *Coef* G' / prod(1:d);



end

