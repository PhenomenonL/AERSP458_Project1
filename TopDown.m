function [U0,U1,U2] = TopDown(alpha,chi)
%TOPDOWN implements top down algorithm to find u and U0 U1 U2
format long;
a0 = 1/2*chi;
b0 = 1;
delta_0 = 1;
tol = 1e-6;
sum_0 = a0/b0;
u0 = a0/b0;
a = zeros(1,200);
b = zeros(1,200);
delta = zeros(1,200);
uArray = zeros(1,200);
sum = zeros(1,200);
a(1) = a0;
b(1) = b0;
delta(1) = delta_0;
uArray(1) = u0;
sum(1) = sum_0;

for n = 2:200
    a(n) = alpha*a0^2; 
    b(n) = b(n-1)+2; %1,3,5,7...
    delta(n) = 1/(1-(a(n)/(b(n-1)*b(n)))*delta(n-1));
    uArray(n) = uArray(n-1)*(delta(n)-1);
    sum(n) = sum(n-1)+uArray(n);
    u = sum(n);    
    if abs(uArray(n)) < tol
        break;
    end
end

U0 = (1-alpha*u^2)/(1+alpha*u^2);
U1 = (2*u)/(1+alpha*u^2);
U2 = (2*u^2)/(1+alpha*u^2);
end

