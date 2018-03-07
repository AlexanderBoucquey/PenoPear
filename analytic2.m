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

syms C1 C2
eqns = [hu*(C1*(exp(-sqrt(V_mu/K_mu)*0.05)/0.05) + C2*(exp(sqrt(V_mu/K_mu)*0.05)/0.05) - C_uamb) == ...
    Dur*(C1*((-sqrt(V_mu/K_mu)*exp(-sqrt(V_mu/K_mu)*0.05)*0.05-exp(-sqrt(V_mu/K_mu)*0.05))/0.05^2) + C2*((sqrt(V_mu/K_mu)*exp(sqrt(V_mu/K_mu)*0.05)*0.05-exp(sqrt(V_mu/K_mu)*0.05))/0.05^2)), hu*(C1*(exp(-sqrt(V_mu/K_mu)*-0.05)/-0.05) + C2*(exp(sqrt(V_mu/K_mu)*-0.05)/-0.05) - C_uamb) == ...
    Dur*(C1*((-sqrt(V_mu/K_mu)*exp(-sqrt(V_mu/K_mu)*-0.05)*-0.05-exp(-sqrt(V_mu/K_mu)*-0.05))/0.05^2) + C2*((sqrt(V_mu/K_mu)*exp(sqrt(V_mu/K_mu)*-0.05)*-0.05-exp(sqrt(V_mu/K_mu)*-0.05))/0.05^2))];
S = solve(eqns, [C1 C2]);
C1 = S.C1;
C2 = S.C2;

r = linspace(-0.05,0.05,500);

Cu = C1*exp(-sqrt(V_mu/K_mu).*r)./r + C2*exp(sqrt(V_mu/K_mu).*r)./r;
plot(r,Cu);