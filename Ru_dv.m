function [ Rudv ] = Ru_dv(  C_u,C_v,V_mu,K_mu,K_mv )
%RU_DV Summary of this function goes here
Rudv = diag(-(V_mu*C_u)./((K_mu+C_u).*(K_mv*(1+C_v/K_mv).^2)));


end

