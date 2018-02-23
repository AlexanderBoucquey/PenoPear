function [ C ] = Cij( r,z )
[area,~,~,~] = triangle(r,z);       
       C= area/60*[[6*r(1)+2*r(2)+2*r(3) 2*r(1)+2*r(2)+r(3) 2*r(1)+r(2)+2*r(3)];
           [2*r(1)+2*r(2)+r(3) 2*r(1)+6*r(2)+2*r(3) r(1)+2*r(2)+2*r(3)];
           [2*r(1)+r(2)+2*r(3) r(1)+2*r(2)+2*r(3) 2*r(1)+2*r(2)+6*r(3)]];
end
