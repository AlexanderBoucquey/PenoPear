function [ Ru ] = R_u( C_u, C_v, V_mu, K_mu, K_mv)
%R_U Summary of this function goes here
Ru = (V_mu.*C_u)./((K_mu + C_u).*(1+C_v./K_mv));


end

