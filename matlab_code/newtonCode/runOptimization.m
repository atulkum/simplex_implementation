%% Demo for optimizing a function using Newton's method
%% we will optimize the Rosenbrock function 100 * (y -x^2)^2 + (1 - x)^2

%% Mesh plot the function
figure(1);
[X,Y]= meshgrid(-15:0.1:15,-15:0.1:15);
Z = 100 * ( Y - X.^2).^2 + (ones(size(X)) - X).^2;
figure(1);
mesh(X,Y,Z);
hold on;

%% Run Newton's method. Use routine rosenbrockDfDDf to get ourselves the gradient and Hessian for the 
%% function to optimize
[x,res,xValues,count] = testNewton(@rosenbrockDfDDf,[5;5], 1000);

%% plot the result on the mesh plot
for i = 1:count
   xValues(i,3) = rosenbrock(xValues(i,:)'); 
end
figure(1);
hold on;
plot3(xValues(:,1), xValues(:,2), xValues(:,3),'--rs')