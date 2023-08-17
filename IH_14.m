clear
clc
r1_vec = [6378 0];
r1 = 6378;
r2 = 42164;
d_theta = linspace(pi/6, pi/4, 1000);
d_time = linspace(7200, 10800, 1000);
r2_vec = [r2*cos(d_theta') r2*sin(d_theta')];
mu = 398600;

c = ((r1 - r2_vec(:,1)).^2 + (r2_vec(:,2)).^2).^.5;
S = .5*(r1 + r2 + c);

% a)
B_m = 2*asin(((S-c)./S).^.5);

d_t_min = (((S.^3)./(mu*8)).^.5) .* (pi - B_m + sin(B_m));
% K=1

a_vals = zeros(1000, 1);
for i = 1:1000
    A = sin(d_theta(i)).*((r1.*r2)./(1-cos(d_theta(i)))).^.5;

    plot_z = linspace(0, 5, 51);
    plot_F = F(r1, r2, d_time(i), d_theta(i), plot_z);
    x_zero = linspace(0, 5, 51);
    F_zero = zeros(1,51);
    %plot(plot_z, plot_F, x_zero, F_zero)

    % Newton's Method
    z = 3.75;
    for j = 1:50
        z = z - (F(r1, r2, d_time(i), d_theta(i), z))/(F_Prime(r1, r2, d_theta(i), z));
    end

    y_new = r1 + r2 + A*((z*stumpS(z) - 1)/(stumpC(z)^.5));
    f = 1-(y_new/r1);
    g = A*(y_new/398600)^.5;
    g_dot = 1-(y_new/r2);

    r1_vec = [r1 0 0];
    r2_vec = [r2*cos(d_theta(i)) r2*sin(d_theta(i)) 0];

    v_1_vec = (1/g).*(r2_vec - f.*r1_vec);

    coe = coe_from_sv(r1_vec, v_1_vec, 398600);

    a_vals(i) = coe(7);
end

tiledlayout(2,2)
nexttile
plot(d_theta, a_vals);
title('Semi-Major Axis vs. Change in True Anomaly')
xlabel('Change in True Anomaly(rad)')
ylabel('Semi-Major-Axis(km)')

% b)

a_m = (r1+r2+c)./4;
nexttile
plot(d_theta, a_m)
title('Minimum Semi-Major Axis vs. Change in  True Anomaly')
xlabel('Change in True Anomaly(rad)')
ylabel('Minimum Semi-Major-Axis(km)')

% c

d_t_parabolic = (2^.5 / (3 * (mu^.5)))*(S.^(3/2) - (sign(sin(d_theta))*(S-c).^(3/2)));
nexttile
plot(d_theta, d_t_parabolic)
title('Parabolic Transfer Time vs. Change in True Anomaly')
xlabel('Change in True Anomaly(rad)')
ylabel('Parabolic Transfer Time (s)')
nexttile
plot(d_theta, d_time)
title('Time to Intercept vs. Change in True Anomaly')
xlabel('Change in True Anomaly(rad)')
ylabel('Time to Intercept (s)')

function F = F(r2, r_2, d_time, d_theta, z)
    A = sin(d_theta).*((r2.*r_2)./(1-cos(d_theta))).^.5;
    F = ((((y(r2, r_2, d_theta, z))./(stumpC(z))).^(3./2)) .* stumpS(z)) + (A.*(y(r2, r_2, d_theta, z)).^.5) - (d_time .* (398600.^.5));
end

function F_Prime = F_Prime(r2, r_2, d_theta, z)
    A = sin(d_theta).*((r2.*r_2)./(1-cos(d_theta))).^.5;
    F_Prime = (1./(2.*(y(r2, r_2, d_theta, z).*(stumpC(z)).^5).^.5)).*(((2.*stumpC(z).*stumpS_Prime(z)) - (3.*stumpC_Prime(z).*stumpS(z))).*y(r2, r_2, d_theta, z).^2 + ((A.*stumpC(z).^(5./2)) + (3.*stumpC(z).*stumpS(z).*y(r2, r_2, d_theta, z))).*y_prime(r2, r_2, d_theta, z));
end

function y = y(r2, r_2, d_theta, z)
    A = sin(d_theta).*((r2.*r_2)./(1-cos(d_theta))).^.5;
    y = r2 + r_2 + (A .* ((z.*stumpS(z) - 1)./(stumpC(z).^.5)));
end

function y_prime = y_prime(r2, r_2, d_theta, z)
    A = sin(d_theta).*((r2.*r_2)./(1-cos(d_theta))).^.5;
    y_prime = (A./4).*(stumpC(z)).^.5;
end

function stumpS_Prime = stumpS_Prime(z)
    stumpS_Prime = (1./(2.*z)).*(stumpC(z) - 3.*stumpS(z));
end

function stumpC_Prime = stumpC_Prime(z)
    stumpC_Prime = (1./(2.*z)).*(1-z.*stumpS(z)-2.*stumpC(z));
end