function [x,y] = solveLPWithInteriorPtMethod(c,A,b)
%% function: solveLPWithInteriorPtMethod
%% Inputs: c,A,b are problem data
%% Problem is assumed to be standard form with
%% max c'*x s.t. Ax <= b, x >=0
%%

[m,n] = size(A);

% Initialize the solutions to some positive values.
x0 = 5*ones(n,1);
y0 = 5*ones(m,1);
xs0 = 5*ones(m,1);
ys0 = 5 * ones(n,1);
% Start with mu = 25
mu=25;
scaleMu=0.1;
M=1E+06;
diff = Inf;
iterCount = 1;
iterLimit=30;
resValues = zeros(100,1);

%% Initialize Primal dual gap.
pdGap = abs(c'*x0 - b' *y0);
tic;
while  (pdGap >= 1E-08 && iterCount < iterLimit)
   %% Convergence criterion is that primal-dual gap must be less than 10^{-15} 
   %% 1. Solve the Log barrier problem for the current value of mu
   [x,y,xs,ys,res] = solveLogBarrier(A,b,c,mu,x0,y0,xs0,ys0);
   %% mu = 1/(m+n)* (x'*ys + y'*xs);
   %% 2. Scale the mu parameter
   mu = scaleMu * mu;
   
   %% record the residual for plotting purposes
   resValues(iterCount,1) = pdGap;
   
   iterCount= iterCount+1;
   
   %% Record the current solution as the previous solution
   x0= x;
   y0 = y;
   xs0 = xs;
   ys0 = ys;
   
   %% Check if norms of x, y and residuals exceed the tolerance M?
   %% If so, diagnose infeasiblitity/unboundess and go away.
   if (norm (x) >= M)
      fprintf ('Warning: Primal seems unbounded and dual infeasible. \n');
      break
   end
   if (norm(y) >= M)
      fprintf('Warning: dual seems unbounded and primal infeasible. \n');
      break
   end

   
   %% Update the primal dual gap.
   
   pdGap = norm(c'*x - b'*y);
end
timeTaken=toc;

figure(1);
plot(1:iterCount-1, resValues(1:iterCount-1));

if (pdGap > 1E-03)
   fprintf ('Primal-Dual Gap did not convege (Final Gap: %f) \n', pdGap);
   fprintf ('Primal Feasibility Gap: %f \n', norm(A*x + xs - b));
   fprintf ('Objectives are: %f PRIMAL, %f DUAL \n', c'*x, b'*y);
   fprintf ('Dual Feasibility Gap: %f \n', norm( A'*y - ys - c));
   fprintf (' KKT-residual: %f  (mu = %f) \n', res,mu);
   return;
end

fprintf (' Optimal solution found with objective value: %f \n', (c'*x));
fprintf (' Number of iterations to converge: %d \n', iterCount);
fprintf ('Primal Feasibility Gap: %f \n', norm(A*x + xs - b));
fprintf ('Dual Feasibility Gap: %f \n', norm( A'*y - ys - c));
fprintf (' Primal-Dual Gap: %f \n', pdGap);
fprintf (' KKT-residual: %f  (mu = %f) \n', res,mu);
fprintf (' Number of Iterations: %d (LIMIT: %d ) \n', iterCount, iterLimit);
fprintf (' Time taken: %f \n', timeTaken);
end