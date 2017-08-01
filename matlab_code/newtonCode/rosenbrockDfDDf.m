function [fx,Jx] = rosenbrockDfDDf(x)
x1 = x(1);
x2 = x(2);

fx=[ -400 * x1* (x2 - x1^2)  + 2 *( x1 -1) ; ...
     2* 100 * (x2 - x1^2) ];
 
dxx = -400 * (x2 - x1)^2 +800 * x1 ^2 + 2;
dxy = -400 * x1;
dyy = 200;

Jx = [ dxx dxy; dxy dyy];
end