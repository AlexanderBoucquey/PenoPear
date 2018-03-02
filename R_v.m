function [ Rv ] = R_v( C_u, C_v, r_q, V_fv, K_mfu,V_mu,K_mu,K_mv )
%R_V Summary of this function goes here
    Rv = r_q.*R_u(C_u,C_v,V_mu,K_mu,K_mv )+V_fv./(1+C_u./K_mfu);


end

