clear variables;
% positieve loop
y = ones(20,1);
x = zeros(20,1);
x(1) = 0.2;
pv = zeros(38,2);
for i=1:19
    y(i+1) = 1-i/18;
    if i < 9
        x(i+1) = x(i) + 0.01;
    elseif i < 14
        x(i+1) = x(i) + 0.03;
    else
        x(i+1) = x(i) - 0.02;
    end
    pv(i,1) = x(i);
    pv(i,2) = y(i);
    pv(39-i,1) = -x(i);
    pv(39-i,2) = y(i);
end
[p,t]=distmesh2d(@dpoly,@huniform,0.1,[-1,-1; 2,1],pv,pv);