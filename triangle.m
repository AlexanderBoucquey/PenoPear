function [area,a,b,c] = triangle(r,z)
    a = circshift(r,2).*circshift(z,1)-circshift(r,1).*circshift(z,2);
    b = circshift(z,2)-circshift(z,1);
    c = circshift(r,1)-circshift(r,2);
    area = r(2)*z(3)+r(1)*z(2) + r(3)*z(1) - r(2)*z(1) - r(1)*z(3) - r(3)*z(2); 
end

