clear
clc
p_0 = [.0107 6748];

[t, h] = ode45(@func, [0, 60*60*24*7], p_0);

plot(t, h(:,1)); % Plotting e vs t
title('Eccentricity Varying With Time')
xlabel('Time(s)')
ylabel('Eccentricity')
%plot(t, h(:,2)); % Plotting a vs t
%title('Semi-Major Axis Varying With Time')
%xlabel('Time(s)')
%ylabel('Semi-Major Axis (Km)')

delta_e = h(1,1) - h(length(h), 1);
delta_a = h(1,2) - h(length(h),2);

function d_dt = func(t, p)
    e = p(1);
    a = p(2);
    
    z_1 = 298;
    density_1 = 1*10^-2;
    z_2 = 442;
    density_2 = 2*10^-3;
    alpha = (log(density_1) - log(density_2))/(z_2 - z_1);
    B = 1*10^-8;
    mu = 398600;
    T = (2*pi*a^(3/2))/(mu^.5);
    M_e = 2*pi*(t/T);
    E = M_e + e*sin(M_e) + .5*e^(2)*sin(2*M_e);
    theta = 2*atan(((1+e)/(1-e))^(.5)*tan(E/2));
    z = a*(1-e*cos(E))- 6378.1;
    density = density_1 * exp(-alpha*(z-z_1));
    F_az = -.5*density*B*(mu/a);
    a_dot = ((2*a^(3/2))/((mu*(1-e^2))^.5))*((1+e*cos(theta))*F_az);
    e_dot = ((a*(1-e^2))/(mu))^.5*(F_az*(cos(theta)*cos(E)));
    d_dt = [e_dot; a_dot];
end