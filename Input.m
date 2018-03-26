function [ ] = Input( mesh_size )
[p,t] = pear_mesh(mesh_size);
rn = boundary(p,0.7);
rp = p(rn,:);

dlmwrite('meshp.txt',p);
dlmwrite('mesht.txt',t);
dlmwrite('rn.txt',rn);
dlmwrite('rp.txt',rp);


end

