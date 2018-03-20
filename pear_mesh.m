function [p,t] = pear_mesh(mesh_size)
y = ones(20,1);
x = zeros(20,1);
k=15;
n=32;
pv = zeros(n,2);
for i=1:n
    %y(i+1) = 1-i/18;
    if i < k
        x(i) = 0.5*sin((i-1)*pi/2/k);
        y(i) = 1.5+0.5*cos((i-1)*pi/2/k);
    elseif i == k
        x(i)=1;
        y(i)=0;
    elseif i <= n
        x(i) = cos((i-k)*pi/2/(n-k));
        y(i) = -sin((i-k)*pi/2/(n-k));
    end
    pv(i,1) = x(i);
    pv(i,2) = y(i);
    
%     pv(n+1-i,1) = -x(i);
%     pv(n+1-i,2) = y(i);
end
pv(n+1,:) = [0 2];
pv = 0.05.*pv;
[p,t]=distmesh2d(@dpoly,@huniform,mesh_size,[-1,-1; 1,2],pv,pv);
end