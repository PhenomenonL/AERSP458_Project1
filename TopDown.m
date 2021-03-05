function [U0,U1,U2] = TopDown(alpha,chi)
%TOPDOWN implements top down algorithm to find u and U0 U1 U2
format long;
a0 = 0.5*chi;
b0 = 1;
delta_0 = 1;

tol = 10e-6; %tolerance
sum0 = a0/b0;
u0 = a0/b0;
a = zeros(1,100);
b = zeros(1,100);
f = zeros(1,100);
delta = zeros(1, 100);
sum = zeros(1,100);
delta(1) = delta_0;
a(1) = a0;
b(1) = b0;
f(1) = u0;
sum(1) = sum0;

for n = 2:100
    if (f(n-1) < tol)
        break;
    end
    a(n) = alpha*a0^2;
    b(n) = b(n-1) + 2;
    delta(n) = 1/(1-(a(n)/(b(n-1)*b(n))*delta(n-1)));
    f(n) = f(n-1)*(delta(n)-1);
    sum(n) = sum(n-1)+f(n);
    u = sum(n);
end

U0 = (1-alpha*u^2)/(1+alpha*u^2);
U1 = (2*u)/(1+alpha*u^2);
U2 = (2*u^2)/(1+alpha*u^2);

end

