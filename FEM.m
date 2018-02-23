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
randn = boundary(p);
randp = p(randn,:);

% Opstellen matrices
for m = 1:length(t)
    r=p(t(m,:),1);
    z=p(t(m,:),2);
    Ku(t(m,:),t(m,:)) = Ku(t(m,:),t(m,:)) + Kij(r,z,Dur,Duz);
    Kv(t(m,:),t(m,:)) = Kv(t(m,:),t(m,:)) + Kij(r,z,Dvr,Dvz);
    C(t(m,:),t(m,:)) = C(t(m,:),t(m,:))+Cij(r,z);
end
l = sqrt((randp(1:end-1,1)-randp(2:end,1)).^2+(randp(1:end-1,2)-randp(2:end,2)).^2);
l(length(randn)) = sqrt((randp(1,1)-randp(end,1)).^2+(randp(1,2)-randp(end,2)).^2);
r = randp(:,1);
K_h = zeros(length(p),length(p));
R_q = zeros(length(p),1);
for o = 1:length(randn)-1
    K_h([randn(o) randn(o+1)],[randn(o) randn(o+1)])= K_h([randn(o) randn(o+1)],[randn(o) randn(o+1)])+Kh(r(o:o+1),l(o));
    R_q([randn(o) randn(o+1)]) = R_q([randn(o) randn(o+1)]) + Rq(r(o:o+1),l(o));
end
  K_h([randn(length(randn)) randn(1)],[randn(length(randn)) randn(1)])= K_h([randn(length(randn)) randn(1)],[randn(length(randn)) randn(1)])+Kh(r([length(randn) 1]),l(length(randn)));
  R_q([randn(length(randn)) randn(1)]) = R_q([randn(end) randn(1)]) + Rq(r([length(randn) 1]),l(length(randn)));

% Lineaire oplossing
C_u0 = (Ku + V_mu/K_mu*C + hu*K_h)\(hu*C_uamb*R_q);
C_v0 = (Kv + hv*K_h)\(r_q*(V_mu/K_mu)*C*C_u0+hv*C_vamb*R_q);
x0 = [C_u0;C_v0];

% Plot lineaire oplossing oxide
figure
[xi,yi] = meshgrid(-1:0.05:1, -1:0.05:2);
zi = griddata(p(:,1),p(:,2),C_u0,xi,yi);
surf(xi,yi,zi);

% Plot lineaire oplossing dioxide
figure
[xi,yi] = meshgrid(-1:0.01:1, -1:0.01:2);
zi = griddata(p(:,1),p(:,2),C_v0,xi,yi);
surf(xi,yi,zi);
% Volledige oplossing
% TODO: Zoek f en fp!
%x = newton(x0,f,fp);