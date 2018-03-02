function [ Rvdu ] = Rv_du( C_u,C_v,r_q,V_fv,K_mfu,V_mu,K_mu,K_mv )
%RV_DU Summary of this function goes here
Rvdu = r_q*Ru_du(C_u,C_v,V_mu,K_mu,K_mv )+ ...
    diag(V_fv*(-1./(K_mfu*(1+C_u/K_mfu).^2)));


end

