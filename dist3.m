function [d] = dist3(x1,y1,z1,x2,y2,z2)
    x=x2-x1;
    y=y2-y1;
    z=z2-z1;
    d=sqrt(x*x+y*y+z*z);
end

