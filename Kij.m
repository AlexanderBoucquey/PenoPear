function [ K ] = Kij( r,z,Dr,Dz )
[area,~,b,c] = triangle(r,z);
K = zeros(3,3);
for i = 1:3
   for j = 1:3       
       K(i,j) = (r(1) + r(2) + r(3))./(12*area).*(Dr*b(i)*b(j) + Dz*c(i)*c(j));
   end
end
end

