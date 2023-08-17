clc 
clear
r_0=6778;
delta_t = 1800;
r_1 = 7378;
delta_theta = 120;
A = sind(delta_theta).*((r_0.*r_1)./(1-cosd(delta_theta))).^.5;


plot_z = linspace(0, 5, 51);
plot_F = F(r_0, r_1, delta_t, delta_theta, plot_z);
x_zero = linspace(0, 5, 51);
F_zero = zeros(1,51);
plot(plot_z, plot_F, x_zero, F_zero)

% Newton's Method
z = 3.75;
for i = 1:50
    z = z - (F(r_0, r_1, delta_t, delta_theta, z))/(F_Prime(r_0, r_1, delta_theta, z));
end


y_new = r_0 + r_1 + A*((z*stumpS(z) - 1)/(stumpC(z)^.5));
f = 1-(y_new/r_0);
g = A*(y_new/398600)^.5;
g_dot = 1-(y_new/r_1);

r_0_vec = [r_0 0 0];
r_1_vec = [r_1*cosd(delta_theta) r_1*sind(delta_theta) 0];

v_0_vec = (1/g).*(r_1_vec - f.*r_0_vec);

coe = coe_from_sv(r_0_vec, v_0_vec, 398600);

r_p = (coe(1)^2/398600)*(1/(1+coe(2)));
h_p = r_p - 6378

function F = F(r_1, r_2, delta_t, delta_theta, z)
    A = sind(delta_theta).*((r_1.*r_2)./(1-cosd(delta_theta))).^.5;
    F = ((((y(r_1, r_2, delta_theta, z))./(stumpC(z))).^(3./2)) .* stumpS(z)) + (A.*(y(r_1, r_2, delta_theta, z)).^.5) - (delta_t .* (398600.^.5));
end

function F_Prime = F_Prime(r_1, r_2, delta_theta, z)
    A = sind(delta_theta).*((r_1.*r_2)./(1-cosd(delta_theta))).^.5;
    F_Prime = (1./(2.*(y(r_1, r_2, delta_theta, z).*(stumpC(z)).^5).^.5)).*(((2.*stumpC(z).*stumpS_Prime(z)) - (3.*stumpC_Prime(z).*stumpS(z))).*y(r_1, r_2, delta_theta, z).^2 + ((A.*stumpC(z).^(5./2)) + (3.*stumpC(z).*stumpS(z).*y(r_1, r_2, delta_theta, z))).*y_prime(r_1, r_2, delta_theta, z));
end

function y = y(r_1, r_2, delta_theta, z)
    A = sind(delta_theta).*((r_1.*r_2)./(1-cosd(delta_theta))).^.5;
    y = r_1 + r_2 + (A .* ((z.*stumpS(z) - 1)./(stumpC(z).^.5)));
end

function y_prime = y_prime(r_1, r_2, delta_theta, z)
    A = sind(delta_theta).*((r_1.*r_2)./(1-cosd(delta_theta))).^.5;
    y_prime = (A./4).*(stumpC(z)).^.5;
end

function stumpS_Prime = stumpS_Prime(z)
    stumpS_Prime = (1./(2.*z)).*(stumpC(z) - 3.*stumpS(z));
end

function stumpC_Prime = stumpC_Prime(z)
    stumpC_Prime = (1./(2.*z)).*(1-z.*stumpS(z)-2.*stumpC(z));
end
