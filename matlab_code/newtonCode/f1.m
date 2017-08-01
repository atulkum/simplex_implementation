function [fx,Jx] = f1(x)

%%
% %%

xx = x(1);
yy = x(2);

z1 = - xx.^2 - 3 * yy.^2 - xx.*yy + 3  * yy + 4 * xx + 5;
z2 = -2 * xx.^2 - 3 * yy.^2 - xx.* yy  + 10 * xx + 3 * yy;

%% the function value
fx=[z1; z2];

%% the Jacobian value
Jx = [ (- 2 * xx - yy + 4)  (-6 *yy - xx + 3 ) ; ...
        (-4 * xx - yy + 10)  (-6 * yy - xx + 3) ];
end