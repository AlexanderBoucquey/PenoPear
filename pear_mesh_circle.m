function [p,t] = pear_mesh_circle(mesh_size)

k=16;
l = 32
n=50;
y = zeros(n+1,1);
x = zeros(n+1,1);
pv = zeros(n+1,2);
for i=2:n+1
    %y(i+1) = 1-i/18;
    if i <= l+1
        x(i) = 0.5*sin((i-1)*pi/2/k);
        y(i) = 0.5*cos((i-1)*pi/2/k);
     elseif i > l+1
        x(i)=0;
        y(i)=-0.5 + (i-l-1)/(n-l);
%     elseif i <= n
%         x(i) = cos((i-k)*pi/2/(n-k));
%         y(i) = -sin((i-k)*pi/2/(n-k));
    end
    pv(i,1) = x(i);
    pv(i,2) = y(i);
    
%     pv(n+1-i,1) = -x(i);
%     pv(n+1-i,2) = y(i);
end
pv(1,:) = [0 0.5];
%pv(n+2,:) = [0 -0.5];
pv = 0.05.*pv
[p,t]=distmesh2d(@dpoly,@huniform,mesh_size,[-1,-1; 1,2],pv,pv);
end