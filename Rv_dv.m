function [ Rvdv ] = Rv_dv( C_u,C_v,r_q,V_mu,K_mu,K_mv )
%RV_DV Summary of this function goes here
Rvdv = r_q.*Ru_dv(C_u,C_v,V_mu,K_mu,K_mv);


end

