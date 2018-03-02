function [ Rudu ] = Ru_du( C_u,C_v,V_mu,K_mu,K_mv )
%RU_DU Summary of this function goes here
Rudu = diag((V_mu*K_mu)./((1+C_v/K_mv).*(K_mu+C_u).^2));

end

