function [ ] = Input( mesh_size, bound_r)
[p,t] = pear_mesh(mesh_size);
rn = boundary(p, bound_r);
rp = p(rn,:);

dlmwrite('meshp.txt',p);
dlmwrite('mesht.txt',t);
dlmwrite('rn.txt',rn);
dlmwrite('rp.txt',rp);


end

