% AERSP 458 Project 1
%% Case1
clear, clc
format long

r0 = [-1.7512, 2.0439, -2.6693]; %[LU]
v0 = [-2.1843, -0.4926, 0.4740]; %[LU/TU]
t0 = 0;
t1 = 1.7; %[TU]
mu = 4*pi^2; %[LU^3/TU^2]

h = cross(r0,v0);
p = norm(h)^2/mu;
a = -mu*norm(r0)/(norm(r0)*norm(v0)^2-2*mu);
e = sqrt(1 - p/a);
alpha = 1/a;
sigma0 = dot(r0,v0)/sqrt(mu);

[chi,U0,U1,U2] = UniversalEqn(alpha, mu, t1, t0, r0, sigma0);

F =  1 - 1/norm(r0)*U2;
G = norm(r0)/sqrt(mu)*U1+sigma0/sqrt(mu)*U2;
r1 = F*r0 + G*v0;
Ft = -sqrt(mu)/(norm(r1)*norm(r0))*U1;
Gt = 1 - U2/norm(r1);
v1 = Ft*r0 + Gt*v0;

energy0 = norm(v0)^2/2 - mu / norm(r0);
energy1 = norm(v1)^2/2 - mu / norm(r1);

h1 = cross(r1,v1);

disp("Case 1: -----------------------------------")
%Part 1: the orbit type
if e == 0
    disp("Orbit shape is CIRCLE");
elseif e == 1
    disp("Orbit shape is PARABOLIC");
elseif e > 1 
    disp("Orbit shape is HYPERBOLIC");
else
    disp("Orbit shape is ELLIPSE");
end
fprintf("Eccentricity e = %f \n", e);
fprintf("Semi-major Axis a = %f \n", a);
fprintf("r1 = \n");
disp(r1');
fprintf("v1 = \n");
disp(v1');
fprintf("Energy at t0 = %f ; Energy at t1 = %f \n", energy0, energy1);
fprintf("h at t0 = \n");
disp(h');
fprintf("h at t1 = \n");
disp(h1');

%% Case2
r0 = [0.6229, 1.3651, -0.0475]; %[LU]
v0 = [8.4232, -1.8123, 4.2091]; %[LU/TU]
t0 = 0;
t1 = 1.2; %[TU]
mu = 4*pi^2; %[LU^3/TU^2]

h = cross(r0,v0);
p = norm(h)^2/mu;
a = -mu*norm(r0)/(norm(r0)*norm(v0)^2-2*mu);
e = sqrt(1 - p/a);
alpha = 1/a;
sigma0 = dot(r0,v0)/sqrt(mu);

[chi,U0,U1,U2] = UniversalEqn(alpha, mu, t1, t0, r0, sigma0);

F =  1 - 1/norm(r0)*U2;
G = norm(r0)/sqrt(mu)*U1+sigma0/sqrt(mu)*U2;
r1 = F*r0 + G*v0;
Ft = -sqrt(mu)/(norm(r1)*norm(r0))*U1;
Gt = 1 - U2/norm(r1);
v1 = Ft*r0 + Gt*v0;

energy0 = norm(v0)^2/2 - mu / norm(r0);
energy1 = norm(v1)^2/2 - mu / norm(r1);

h1 = cross(r1,v1);

disp("Case 2: -----------------------------------")
%Part 1: the orbit type
if e == 0
    disp("Orbit shape is CIRCLE");
elseif e == 1
    disp("Orbit shape is PARABOLIC");
elseif e > 1 
    disp("Orbit shape is HYPERBOLIC");
else
    disp("Orbit shape is ELLIPSE");
end
fprintf("Eccentricity e = %f \n", e);
fprintf("Semi-major Axis a = %f \n", a);
fprintf("r1 = \n");
disp(r1');
fprintf("v1 = \n");
disp(v1');
fprintf("Energy at t0 = %f ; Energy at t1 = %f \n", energy0, energy1);
fprintf("h at t0 = \n");
disp(h');
fprintf("h at t1 = \n");
disp(h1');

