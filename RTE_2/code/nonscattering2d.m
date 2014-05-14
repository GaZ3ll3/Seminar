function x=  TEST( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global data

N = 41;
Sigma = zeros(N,N);

% sigma = 0.1 + x^3 + y^3


for i = 1:N
    for j = 1: N
        Sigma(i,j) = sin(pi*(j-1)/(N-1))*sin(pi*(i-1)/(N-1));
    end
end


data = zeros(N,N);


for i = 1:N
    for j = 1: N
        data(i,j) = (sum(Sigma(i,1:j)) - sum(Sigma(1:i,j)));
    end
end


options = optimoptions('fminunc','GradObj','on','DerivativeCheck','off', 'Algorithm','quasi-newton','MaxFunEvals',80000,'Display','iter','TolFun',1e-10,'TolX',1e-10,'MaxIter',1000);


x0 = zeros(N^2,1);
x = fminunc(@opt,x0,options);

figure(1);
surf(Sigma - reshape(x0,N,N));
figure(2);
surf(Sigma - reshape(x,N,N));


d = reshape(x,N,N);
mydata = zeros(N,N);

for i = 1:N
    for j = 1:N
        mydata(i,j) = (sum(d(i,1:j)) - sum(d(1:i,j)));
    end
end

end


function [ret,g] = opt(sigma)
global data


n = size(sigma,1);
N = sqrt(n);
d = reshape(sigma,N,N);

ret = 0.0;

mydata = zeros(N,N);
grad = zeros(N,N);

for i = 1:N
    for j = 1:N
        mydata(i,j) = (sum(d(i,1:j)) - sum(d(1:i,j)));
    end
end


for i =1:N
    for j = 1:N
        ret = ret + (data(i,j) - mydata(i,j))^2;
    end
end

for i = 1:N
    for j = 1:N
        grad(i,j) = 2*sum(data(i:end,j) - mydata(i:end,j)) - 2*(sum(data(i,j:end) - mydata(i,j:end)));
    end
end

grad(1,1:N) = 0;
grad(N,1:N) = 0;
grad(1:N,N) = 0;
grad(1:N,1) = 0;
g = reshape(grad,N^2,1);


end
