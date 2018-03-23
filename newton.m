function [Cu,Cv,N] = newton(C_u0,C_v0,F,J)
% newton-raphson algorithm
N =1000; eps = 1.e-6; % define max. no. iterations and error

u = size(C_u0,1);
Cu = C_u0;
Cv = C_v0;
x = [Cu;Cv];
res = eps + 1;

while (N>0)
 N = N - 1;
 disp(res)
 disp(N)
 xn = x - J(Cu,Cv)\(F(Cu,Cv));
 res = norm(x-xn);
if res < eps
 disp('converged');
 Cu = xn(1:u);
 Cv = xn(u+1:end);
 return;
end 
 x = xn;
 Cu = x(1:u);
 Cv = x(u+1:end);
end
error('No convergence');
end
% end function 

