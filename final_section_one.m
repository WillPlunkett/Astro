clear 
clc

% Initialize constants and setup params
mu = 398600;

r1 = 200+6378;
theta_A_0 = 0;
r3 = 4000 + 6378;
theta_B_0 = pi/2;

% Begin by trying Hohmann transfer equations.
% Calculate energy of the three orbits
h1 = (mu*r1)^.5;
h2 = ((2*mu)^.5)*((r1*r3)/(r1+r3))^.5;
h3 = (mu*r3)^.5;

e2 = 1-((h2^2)/(mu*r3));

% Calcualte required speeds at intersection points
v_A_1 = h1/r1;
v_A_2 = h2/r1;
delta_v1 = v_A_2 - v_A_1;

v_B_2 = h2/r3;
v_B_3 = h3/r3; 
delta_v2 = v_B_3 - v_B_2;

% Total delta v required is less than specified amount.
delta_v_tot = delta_v1 + delta_v2;

a2 = (r1+r3)/2;
T2 = ((2*pi)/(mu^.5))*a2^(3/2);
% Transfer will take half the ellipse period time
transfer_time = T2/2;
% If we want the entire ordeal to take <8000s, the latest t_1 is below.
latest_launch = 8000 - transfer_time;

% Hohmann transfer meets requirements, now we need to find out if a launch
% time within this window exists that will allow the two objects to
% intersect

T1 = ((2*pi)/(mu^.5))*r1^(3/2);
T3 = ((2*pi)/(mu^.5))*r3^(3/2);

% Calculate satellite location (angles) for 1000 times between 0 and latest
% time determined above
t_wait = linspace(0, latest_launch, 1000);
A_angles = rad_angle_A(t_wait, T1);
B_angles = rad_angle_B(t_wait + transfer_time, T3);

% Graph
%plot(t_wait, A_angles, t_wait, B_angles)
%legend('A angles','B Angles')
%xlabel('time (s)')
%ylabel('Angle (rad)')
%title('Angular Posision of Spacecraft as a Function of Time')

% Inspecting, t ~ 1250 seems to be where they collide.

% Solving analytically, we get the following for the wait time
t_analytic = ((2*pi*transfer_time/T3)-(pi/2))/((2*pi/T1)-(2*pi/T3));
t_intercept = t_analytic + transfer_time;

diff = rad_angle_A(t_analytic, T1) - rad_angle_B(t_analytic + transfer_time, T3);

% Testing the solution numerically.
rA =[r1*cos(2*pi*t_analytic/T1) r1*sin(2*pi*t_analytic/T1) 0];
r_dot_A = [-v_A_2*sin(2*pi*t_analytic/T1) v_A_2*cos(2*pi*t_analytic/T1) 0];
fA = [rA r_dot_A];

[t_num_A, output_A_num] = ode45(@func, [0, transfer_time], fA);
r_final = (output_A_num(53,1)^2 + output_A_num(53,2)^2)^.5;

rB =[r3*cos((pi/2) + 2*pi*t_analytic/T3) r3*sin((pi/2) + 2*pi*t_analytic/T3) 0];
r_dot_B = [-v_B_3*sin((pi/2) + 2*pi*t_analytic/T3) v_B_3*cos((pi/2) + 2*pi*t_analytic/T3) 0];
fB = [rB r_dot_B];

[t_num_B, output_B_num] = ode45(@func, [0, transfer_time], fB);

% Graphing earth to prove it doesn't hit it
earth_x = linspace(-6378, 6378, 10000);
earth_y_top = (6378.^2 - earth_x.^2).^.5;
earth_y_bottom = -earth_y_top;
earth_x_tot = cat(2, earth_x, earth_x);
earth_y_tot = cat(2, earth_y_top, earth_y_bottom);

plot(output_A_num(:, 1), output_A_num(:, 2), output_B_num(:, 1), output_B_num(:, 2), earth_x_tot, earth_y_tot);
legend('A Orbit', 'B Orbit', 'Earth')
title('Orbit of Crafts A and B Between First Delta V and Intersection')


disp('Total delta v:')
disp(delta_v_tot)
disp('Time of first rocket fire:')
disp(t_analytic)
disp('Time of rendezvous:')
disp(t_intercept)
disp('Eccentricity of transfer:')
disp(e2)
disp('Angular Momentum of transfer:')
disp(h2)
disp('r_p of transfer:')
disp(r1)


% Orbital EOM
function drdt = func(t, f)
    r_x = f(1);
    r_y = f(2);
    r_z = f(3);
    v_x = f(4);
    v_y = f(5);
    v_z = f(6);
    r = norm([r_x, r_y, r_z]);
    a_x = ((-398601.9.*r_x)./(r.^3));
    a_y = ((-398601.9.*r_y)./(r.^3));
    a_z = ((-398601.9.*r_z)./(r.^3));
drdt = [v_x ; v_y ; v_z ; a_x ; a_y ; a_z];
end


% Gets location of satellite A at time t
function val = rad_angle_A(t_wait, T1)
    val = ((2.*pi.*t_wait)./T1) + pi;
end

% Gets location of satellite B at time t
function val = rad_angle_B(t, T3)
    val = (pi/2) + ((2.*pi.*t)./T3);
end