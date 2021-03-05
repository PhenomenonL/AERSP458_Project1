% AERSP 458 Project 1
%% Case1
clear, clc

r0 = [-1.7512, 2.0439, -2.6693]; %[LU]
v0 = [-2.1843, -0.4926, 0.4740]; %[LU/TU]
t0 = 0;
t1 = 1.7; %[TU]
mu = 4*pi^2; %[LU^3/TU^2]

h = cross(r0,v0);
p = norm(h)^2/mu;
a = -mu*norm(r0)/(norm(r0)*norm(v0)^2-2*mu);
e = sqrt(1 - p/a);

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




%% Case2
clear,clc
r0 = [0.6229, 1.3651, -0.0475]; %[LU]
v0 = [8.4232, -1.8123, 4.2091]; %[LU/TU]
t0 = 0;
t1 = 1.2; %[TU]
mu = 4*pi^2; %[LU^3/TU^2]

h = cross(r0,v0);
p = norm(h)^2/mu;
a = -mu*norm(r0)/(norm(r0)*norm(v0)^2-2*mu);
e = sqrt(1 - p/a);

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

