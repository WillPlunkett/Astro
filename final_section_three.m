clear
clc
mu1 = 5959.9;
mu_j = 126686531.9;
% radius of orbit, not of moon
r_io = 421800;
h1=50;
delta_v = 2.4719;
r1 = 1871.5;
% v1 is the ship velocity in parking orbit
v1 = (mu1/r1)^.5;
% v2 is ship velocity after impulsive thrust
v2 = v1 + delta_v;
% v_esc is the escape velocity from IO
v_esc = (2*mu1/r1)^.5;
% v_inf is the velocity the craft escapes IO's SOI with
v_inf = (v2^2 - v_esc^2)^.5;

% vio is the velocity of IO around Jupiter
v_io = (mu_j/r_io)^.5;
% v3 transitions to Jupiter as the dominating gravitational body, adding
% the velocity of IO and the craft together.
v3 = v_io + v_inf;

% Hohmann transfer problem now with r_p = r_io, and v_0 = v3
a = mu_j / -(2*(((v3^2)/2)-(mu_j/r_io)));
r_a = (2*a) - r_io;
% This corresponds with the orbital radius of Ganymede.

% Orbital radius, mu, and velocity of Ganymede
mu_g = 9887.8;
r_g = 1070400;
v_g = (mu_j/r_g)^.5;

% velocity of craft at apoapsis of elliptical orbit around Jupiter
v4 = (2*(-(mu_j/(2*a))+(mu_j/r_a)))^.5;

% Since we can ony slow down by 2 km/s, this leaves ~.7 km/s of velocity
delta_v_slow = v_g - v4;

v_capture = delta_v_slow - 2;

% This is a very large orbit, but they don't have to expell all of their
% fuel, and decreasing the final impulsive thrust would allow a tighter
% orbit (Some thrust is needed, however, as without it, the orbit is below
% the surface of the moon.
r_capture = mu_g/(v_capture^2);
