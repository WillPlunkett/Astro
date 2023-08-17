
a = 7016;
e = .05;
RA = 0;
incl = 45 * pi / 180;
w = 20 * pi / 180;
TA = 10 * pi / 180;
mu = 398600;
h = (a*(1-e^2)*mu)^.5;

coe = [h e RA incl w TA];
[r, v] = sv_from_coe(coe, mu);