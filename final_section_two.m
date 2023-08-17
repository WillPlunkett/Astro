clear
clc

% Define constants, and calculate orbital params
r_p = 3396 + 60;
r_a = 3396 + 5000;
a = (r_p + r_a)./2;
e = 1-(r_p./a);
R = 3396;
mu = 42828;
J_2 = 1.95545 .* (10.^-3);
J_4 = -1.53774 .* (10.^-5);
n = (mu.^.5)./(a.^(3./2));
p = a.*(1-e.^2);
i = asind((4./5).^.5);


omega_dot_1 = omega_dot_one(i);
omega_dot_2 = omega_dot_two(i);
omega_dot = omega_dot_1 + omega_dot_2;
omega_change_yearly = omega_dot .* 60 .* 60 .* 24 .* 687.05;

% Find i that results in constant omega_dot by solving equations. Or,
% approx.

% Plotting to see i-omega_dot curve. Intercept appears at i ~ 63
test_i = linspace(0, 90, 1000);
test_omega_dot = omega_dot_one(test_i) + omega_dot_two(test_i);
plot(test_i, test_omega_dot)
xlabel('i (deg)');
ylabel('Omega dot (degrees/s)')
title('Change in Argument of Perigee With Respect To Inclination')

% Use bisection method to find the root
i_final = bisect(@omega_dot_tot, 62, 64);


function val = omega_dot_tot(i)
    val = omega_dot_one(i) + omega_dot_two(i);
end

function val = omega_dot_one(i)
    r_p = 3396 + 60;
    r_a = 3396 + 5000;
    a = (r_p + r_a)./2;
    e = 1-(r_p./a);
    R = 3396;
    mu = 42828;
    J_2 = 1.95545 .* (10.^-3);
    n = (mu.^.5)./(a.^(3./2));
    p = a.*(1-e.^2);
    val = -((3.*n.*J_2.*R.^2)./(2.*p.^2)).*((5./2).*sind(i).^2 - 2);
end

function val = omega_dot_two(i)
    r_p = 3396 + 60;
    r_a = 3396 + 5000;
    a = (r_p + r_a)./2;
    e = 1-(r_p./a);
    R = 3396;
    mu = 42828;
    J_2 = 1.95545 .* (10.^-3);
    J_4 = -1.53774 .* (10.^-5);
    n = (mu.^.5)./(a.^(3./2));
    p = a.*(1-e.^2);
    val = ((9.*n.*(J_2.^2).*R.^4)./(p.^4)).*((4+((7./12).*e.^2)+(2.*(1-e.^2).^.5))-(sind(i).^2 .* ((103./12)+((3./8).*e.^2)+((11./2).*(1-e.^2).^.5)))+((sind(i).^4).*((215./48)-((15./32).*e.^2)+((15./4).*(1-e.^2).^.5)))-((35.*J_4)./(18.*J_2.^2)).*(((12./7)+(27./14).*e.^2)-((sind(i).^2).*((93./14)+(27./4).*e.^2))+((sind(i).^4).*((21./4)+(81./16).*e.^2))));
end