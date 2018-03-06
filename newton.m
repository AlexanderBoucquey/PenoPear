function [Cu,Cv,N] = newton(C_u0,C_v0,F,J)
% newton-raphson algorithm
N = 5; eps = 1.e-12; % define max. no. iterations and error
maxval = 10000.0; % define value for divergence
u = size(C_u0);
v = size(C_v0);
n = norm(F(C_u0,C_v0));
xx = [C_u0;C_v0];
Cu = xx(1:u);
Cv = xx(u+1:u+v);
while (N>0)
 dx = J(Cu,Cv)\(-F(Cu,Cv));
 xn = xx + dx;
% xn = xx-F(xx)/J(xx);
if (norm(F(xn(1:u),xn(u+1:u+v)))/n)<eps
 disp('converged');
 x=xn;iter=100-N;
 Cu = x(1:u);
 Cv = x(u+1:u+v);
 return;
end
if norm(F(Cu,Cv))>maxval
 disp(['iterations = ',num2str(iter)]);
 error('Solution diverges');
 break;
end
    
 N = N - 1
 xx = xn;
 Cu = xx(1:u);
 Cv = xx(u+1:u+v);
end
error('No convergence');
end
% end function 

