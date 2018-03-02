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
V_mu = 2.39E10-4*exp(80200/R_g*(1/293.15-1/T));
V_fv = 1.61E10-4*exp(56700/R_g*(1/293.15-1/T));
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

% Mesh
[p,t] = pear_mesh();
Ku = zeros(length(p),length(p));
Kv = zeros(length(p),length(p));
C = zeros(length(p),length(p));

% Randen
rn = boundary(p);
rp = p(rn,:);

% Opstellen matrices
for m = 1:length(t)
    r=p(t(m,:),1);
    z=p(t(m,:),2);
    Ku(t(m,:),t(m,:)) = Ku(t(m,:),t(m,:)) + Kij(r,z,Dur,Duz);
    Kv(t(m,:),t(m,:)) = Kv(t(m,:),t(m,:)) + Kij(r,z,Dvr,Dvz);
    C(t(m,:),t(m,:)) = C(t(m,:),t(m,:)) + Cij(r,z);
end
l = sqrt((rp(1:end-1,1)-rp(2:end,1)).^2+(rp(1:end-1,2)-rp(2:end,2)).^2);
l(length(rn)) = sqrt((rp(1,1)-rp(end,1)).^2+(rp(1,2)-rp(end,2)).^2);
r = rp(:,1);
K_h = zeros(length(p),length(p));
R_q = zeros(length(p),1);

for o = 1:length(rn)-1
    if((r(o)>1E-13) || (r(o+1)>1E-13)) 
        K_h([rn(o) rn(o+1)],[rn(o) rn(o+1)])= K_h([rn(o) rn(o+1)],[rn(o) rn(o+1)])+Kh(r(o:o+1),l(o));
       % disp(o);
        R_q([rn(o) rn(o+1)]) = R_q([rn(o) rn(o+1)]) + Rq(r(o:o+1),l(o));
    end
end


if((r(end)>1E-13) || (r(1)>1E-13)) 
  K_h([rn(end) rn(1)],[rn(length(rn)) rn(1)])= K_h([rn(end) rn(1)],[rn(end) rn(1)])+Kh(r([length(rn) 1]),l(length(rn)));
  R_q([rn(end) rn(1)]) = R_q([rn(end) rn(1)]) + Rq(r([length(rn) 1]),l(length(rn)));
end


%Lineaire oplossing
C_u0 = (Ku + V_mu/K_mu*C + hu*K_h)\(hu*C_uamb*R_q);
C_v0 = (Kv + hv*K_h)\(r_q*(V_mu/K_mu)*C*C_u0+hv*C_vamb*R_q);

x0 = [C_u0;C_v0];

% Plot lineaire oplossing oxide
figure
[xi,yi] = meshgrid(-0.05:0.001:0.05, -0.05:0.001:0.1);
zi = griddata(p(:,1),p(:,2),C_u0,xi,yi);
surf(xi,yi,zi);

% Plot lineaire oplossing dioxide
figure
[xi,yi] = meshgrid(-0.05:0.001:0.05, -0.05:0.001:0.1);
zi = griddata(p(:,1),p(:,2),C_v0,xi,yi);
surf(xi,yi,zi);

figure
scatter(rp(:,1),rp(:,2));

% Volledige oplossing
% TODO: Zoek f en fp!
F = @(C_u,C_v) [Ku*C_u+C*Ru(C_u,C_v,V_mu,K_mu,K_mv)+hu*(K_h*C_u-R_q*C_uamb);...
    Kv*C_v-C*Rv(C_u,C_v,r_q,V_fv,K_mfu,V_mu,K_mu,K_mv)+hv*(K_h*C_v-R_q*C_vamb)];
J = @(C_u,C_v) [[Ku+C*Rudu(C_u,C_v, V_mu,K_mu,K_mv)+hu*K_h...
    C*Rudv(C_u,C_v, V_mu,K_mu,K_mv)];...
    [-C*Rvdu(C_u,C_v,r_q,V_mfv,K_mfu,V_mu,K_mu,K_mv)...
    Kv-C*Rvdv(C_u,C_v,r_q,V_mu,K_mu,K_mv)+hv*K_h]];


[C_u,C_v] = newton(x0,F,J);