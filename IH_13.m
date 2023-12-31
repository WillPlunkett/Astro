clear 
clc
Omega_dot_ss = 7.292115*10^-5 * pi/180;
J2 = 1.08263*10^-3;
e = [0 0 0 0 0 0 0];
N = [11 12 13 14 15 16 17];
mu = 398600;
R = 6378;
a = (86400 * mu^.5 ./ (N.*pi*2)).^(2/3);
%i = acosd((a.^(7/2).*Omega_dot_ss.*2.*(1-e.^2).^2)/ (-3.*J2.*(mu.^.5).*R.^2));
T = (2.*pi.*a.^(3/2))./ (mu).^.5;
track_len_deg = T.*360 / (24*60*60);
track_len_km = track_len_deg.*111;
h = (a.*mu.*(1-e.^2)).^.5;
r_a = (h.^2)./(mu*(1-e));

N_2 = [24 25 26 27 28 29 30 31 32];
revs_day = N_2./2;
a_2 = (86400 * mu^.5 ./ (revs_day.*pi*2)).^(2/3);
T_2 = (2.*pi.*a_2.^(3/2))./ (mu).^.5;
track_len_deg_2 = T_2.*360 / (24*60*60);
track_len_km_2 = track_len_deg_2.*111;

N_3 = [32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50];
revs_day_3 = N_3./3;
a_3 = (86400 * mu^.5 ./ (revs_day_3.*pi*2)).^(2/3);
T_3 = (2.*pi.*a_3.^(3/2))./ (mu).^.5;
track_len_deg_3 = T_3.*360 / (24*60*60);
track_len_km_3 = track_len_deg_3.*111;


