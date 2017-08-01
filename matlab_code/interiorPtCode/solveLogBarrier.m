function [x,y,xs,ys, res]= solveLogBarrier(A,b,c,mu,x0,y0,xs0, ys0)
%% Function: solveLogBarrier
%% A, b, c: Problem Data
%% We are solving a standard form problem.
%% max c' x s.t. A x <= b and x >= 0
%% 
%% x0, y0, xs0, ys0: Initial values. Must be > 0 for log barrier problem.
%% mu: Value of mu for which we solve.


%% Initialize values
[m,n] = size(A);
x = x0;
y = y0;
xs = xs0;
ys = ys0;
MAXCOUNT=50; %% MAX number of interations
resTOL=1E-08; %% desired tolerance.
res = Inf;
count = 1;
rho = 0.8;

while ( count < MAXCOUNT && res >= resTOL) %% Main loop for Newton Iteration
    %% Compute F(Xi)
    F = [ A * x + xs - b; % Primal infeasibility
        A' * y - ys - c; % dual infeasibility
        diag(x) * diag(ys) * ones(n,1) - mu*ones(n,1); % mu complementary gap
        diag(y) * diag(xs)*ones(m,1) - mu *ones(m,1)];
    %% Compute its derivative.
    dF = [ A     zeros(m,m) eye(m) zeros(m,n) ;
        zeros(n,n) A' zeros(n,m) -eye(n);
        diag(ys) zeros(n,m) zeros(n,m) diag(x) ;
        zeros(m,n) diag(xs) diag(y) zeros(m,n)];
    
    if (abs(rcond(dF)) < 1E-16)
       fprintf('Condition number of Hessian is not very good.PRIMAL infeasible or DUAL infeasible. \n' ); 
       %%return;
    end
    %% The Newton Step: Solving 
    delta = -dF\F;
    
    %% Compute the scale factor lambda
    lambdaInv = max( abs(delta) ./ [x; y; xs; ys] );
    lambda = rho /lambdaInv;
    
    %% Make sure that we are not fooled into taking a 
    %% large step in the delta direction.
     if (lambda >= 1)
         lambda = 1;
     end
    count = count +1;
    
    %% Make the Update to Current values 
    x = x + lambda * delta(1:n,1);
    y = y + lambda * delta(n+1:n+m,1);
    xs = xs + lambda * delta(n+m+1:n+2*m,1);
    ys = ys + lambda * delta(n+2*m+1:2*n+2*m,1);
    
    %% Is everything still OK??
    assert(min(x) > 0);
    assert(min(xs)> 0);
    assert(min(y) > 0);
    assert(min(ys) > 0);
%      mu = 1/(m+n) * (sum(diag(y) *xs)+ sum(diag(x)*ys));     
%      mu = rho *mu;
    res = norm(F);
end



end