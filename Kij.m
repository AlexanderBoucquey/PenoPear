function [ K ] = Kij( p,t )
[area,~,b,c] = triangle(r,z);
for i = 1:length(t,1)
   for j = 1:length(p,1)
       kij = 0;
       K(i,j) = 0;
   end
end
end

