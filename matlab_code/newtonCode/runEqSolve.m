%% Plotting commands for generating plot in slides.
figure(1);
[X,Y]= meshgrid(-15:0.1:15,-15:0.1:15);
Z1= - X.^2 - 3 * Y.^2 - X.*Y + 3  * Y + 4 * X + 5;
Z2 = -2 * X.^2 - 3 * Y.^2 - X.* Y  + 10 * X + 3 * Y;
contour(X,Y,Z1,'ShowText','on');
hold on;
contour(X,Y,Z2,'ShowText','on');
hold on;
[x,res] = testNewton(@f1,[20;-20], 1000)