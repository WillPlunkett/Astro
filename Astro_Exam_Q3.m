f_0 = [-.01 .001 .02 0 .00003 0];
I_1 = 40;
I_2 = 46;
I_3 = 30;
h_2 = .3;
w_0 = .00105;
[t, h] = ode45(@func, [0, 60*70], f_0);
plot(t, h(:, 1), t, h(:, 2), t, h(:, 3))
title('Yaw, Pitch, and Roll as a Function of Time')
xlabel('Time (s)')
ylabel('Angle (rads)')
legend('Yaw', 'Pitch', 'Roll')
function d_dt = func(t,f)
yaw = f(1);
pitch = f(2);
roll = f(3);
yaw_dot = f(4);
pitch_dot = f(5);
roll_dot = f(6);
I_1 = 40;
I_2 = 46;
I_3 = 30;
h_2 = .3;
w_0 = .00105;
yaw_dot_dot = ((-(((w_0^2)*(I_2-I_1))+(w_0*h_2))*yaw)+(((w_0*(I_1-I_2+I_3))-h_2)*roll_dot))/I_3;
pitch_dot_dot = 0;
roll_dot_dot = ((-(((w_0^2)*(I_2-I_3))+(w_0*h_2))*roll)-(((w_0*(I_1-I_2+I_3))-h_2)*yaw_dot))/I_1;
d_dt = [yaw_dot; pitch_dot; roll_dot; yaw_dot_dot; pitch_dot_dot; roll_dot_dot];
end