clear all
close all
Dur = 2.8E-10;
Duz = 1.1E-9;
Dvr = 2.32E-9;
Dvz = 6.97E-9;
R_g = 8.314;
T = 273;
syms V_mu(T) V_mfv(T) C_uamb(T,nuu) C_vamb(T,nuv)
V_mu(T) = 2.39E10-4*exp(80200/R_g*(1/293.15-1/T));
V_fv(T) = 1.61E10-4*exp(56700/R_g*(1/293.15-1/T));
K_mu = 0.4103;
K_mv = 27.2438;
K_mfu = 0.1149;
r_q = 0.97;
hu = 7E-7;
hv = 7.5E-7;
p_atm = 101300;
C_uamb(T,nuu) = p_atm*nuu/(R_g*T);
C_vamb(T,nuv) = p_atm*nuv/(R_g*T);
[p,t] = pear_mesh();
Ku = zeros(length(p),length(p));
Kv = zeros(length(p),length(p));
C = zeros(length(p),length(p));

for m = 1:length(t)
    r=p(t(m,:),1);
    z=p(t(m,:),2);
    Ku(t(m,:),t(m,:)) = Ku(t(m,:),t(m,:)) + Kij(r,z,Dur,Duz);
    Kv(t(m,:),t(m,:)) = Kv(t(m,:),t(m,:)) + Kij(r,z,Dvr,Dvz);
    C(t(m,:),t(m,:)) = C(t(m,:),t(m,:))+Cij(r,z);
end
% Lineaire oplossing
% TODO: x0 = lineaire oplossing
x0 = 0;

% Volledige oplossing
% TODO: Zoek f en fp!
%x = newton(x0,f,fp);