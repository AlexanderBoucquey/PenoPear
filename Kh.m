function [ K_h ] = Kh( r,l )       
       K_h= (l/12).*[[3*r(1)+r(2) r(1)+r(2)];
            [r(1)+r(2) r(1)+3*r(2)]];
end