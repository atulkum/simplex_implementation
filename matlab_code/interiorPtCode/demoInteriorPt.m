m = 100;
n = 110;
A = randn(m,n);
b = 10*A*abs(randn(n,1));
c = randn(n,1);

%%tic;
[x,y] = solveLPWithInteriorPtMethod(c,A,b);
%%toc;

fprintf ('Trying linprog...\n');
tic;
[xP,fVal] = linprog(-c,A,b,[],[], zeros(n,1));
tElapsed = toc;

fprintf('Linprog gave optimum: %f with time taken %f \n', -fVal, tElapsed);

