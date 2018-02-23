function [p,t] = pear_mesh()
y = ones(20,1);
x = zeros(20,1);
k=13;
n=54;
pv = zeros(n,2);
for i=1:n/2
    %y(i+1) = 1-i/18;
    if i < k
        x(i) = 0.5*sin((i-1)*pi/2/k);
        y(i) = 1.5+0.5*cos((i-1)*pi/2/k);
    elseif i == k
        x(i)=1;
        y(i)=0;
    elseif i <= n/2
        x(i) = cos((i-k)*pi/2/(n/2-k-1));
        y(i) = -sin((i-k)*pi/2/(n/2-k-1));
    end
    pv(i,1) = x(i);
    pv(i,2) = y(i);
    pv(n+1-i,1) = -x(i);
    pv(n+1-i,2) = y(i);
end

[p,t]=distmesh2d(@dpoly,@huniform,0.001,[-1,-1; 1,2],pv,pv);
end