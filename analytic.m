clear all
close all

% Constanten
Dur = 2.8E-10;
Duz = 1.1E-9;
Dvr = 2.32E-9;
Dvz = 6.97E-9;
R_g = 8.314;
T = 272.15;
nuu = 2/100;
nuv =0.7/100;
V_mu = 2.39E-4*exp(80200/R_g*(1/293.15-1/T));
V_mfv = 1.61E-4*exp(56700/R_g*(1/293.15-1/T));
K_mu = 0.4103;
K_mv = 27.2438;
K_mfu = 0.1149;
r_q = 0.97;
hu = 7E-7;
hv = 7.5E-7;
p_atm = 101300;

% Uitrekenbare constanten
C_uamb = p_atm*nuu/(R_g*T);
C_vamb = p_atm*nuv/(R_g*T);

%syms u(r) v(r) z(r) w(r)

% ode1 = diff(u) == z./r.^2;
% ode2 = diff(v) == w./r.^2;
% ode3 = diff(z) == r.^2.*(V_mu.*u./((K_mu+u).*(1+v./K_mv)));
% ode4 = diff(w) == -r.^2.*r_q.*(V_mu.*u./((K_mu+u).*(1+v./K_mv)))-r.^2.*(V_mfv./(1+u./K_mfu));
% odes = [ode1; ode2; ode3; ode4];
[r,y] = ode45(@(r,y)[y(3)./r.^2;y(4)./r.^2;r.^2.*(V_mu.*y(1)./((K_mu+y(1)).*(1+y(2)./K_mv)));-r.^2.*r_q.*(V_mu.*y(1)./((K_mu+y(1)).*(1+y(2)./K_mv)))-r.^2.*(V_mfv./(1+y(1)./K_mfu));],[0.1 0],[0.1;0.6;0;0]);

% S = dsolve(odes);
% 
% uSol(r) = S.u;
% vSol(r) = S.v;
% zSol(r) = S.z;
% wSol(r) = S.w;
plot(r,y(:,1),r,y(:,2));
