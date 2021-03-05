function [chi,U0,U1,U2] = UniversalEqn(alpha, mu, t, t0, r0, sigma0)
%UNIVERSALEQN iterativly solve the Universal Time Eqn
tol = 1e-6;
chiOld = alpha*sqrt(mu)*(t-t0); %initial chi value Guess
[U0,U1,U2] = TopDown(alpha, chiOld); %inital U0, U1, U2
chiDiff = abs(chiOld);

while chiDiff >= tol
    fX = chiOld - sigma0 - sigma0*U0 - (1-alpha*r0)*U1-alpha*sqrt(mu)*(t-t0); %f(X_old)
    dfX = 1+sigma0*alpha*U1 - (1-alpha*r0)*U0; %f'(X_old)
    chiNew = chiOld - fX/dfX; %Newton-Raphson
    chiDiff = abs(chiNew-chiOld); % error
    chiOld = chiNew; %re-intialize chiOld
    [U0,U1,U2] = TopDown(alpha, chiNew); %calculate U0,U1,U2 with chiNew
end
 chi = chiNew;
   
end

