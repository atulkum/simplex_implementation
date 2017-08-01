function [x,res,xValues,count] = testNewton(f,x0,maxIters)


count = 1;
x = x0;
res=Inf;
[n, ~] = size(x);
xValues=zeros(maxIters+1,n);
resValues=zeros(maxIters+1,1);
while (count <= maxIters && res >= 1E-10)
    %% Evaluate function and its derivative
    [fx,J] = f(x);
    %% Compute Newton step.
    delta = -J\fx;
   
    
    
    %% Update count and residue
    res = norm(fx);
    %% Record keeping for plotting. 
    %% Can be commented out.
    xValues(count,:) = x';
    
    resValues(count,1) = res;

     %% update current value of x
    x = x + delta;
    count=count+1;
end


xValues(count,:) = x';
resValues(count,1) = res;

if (res >= 1E-10)
   fprintf (' Warning: did not converge to tolerance in %d iterations -- achieved residue is %f \n', maxIters, res); 
else
    fprintf (' Newton converged to residue: %f in %d steps \n', res, count);
    figure(1);
    hold on;
    plot3(xValues(1:count,1), xValues(1:count,2),zeros(count,1),'--rs','LineWidth',2,...
                        'MarkerEdgeColor','k',...
                        'MarkerFaceColor','g')
    xlabel('x1');
    ylabel('x2');
    title('Newton method iterates')
    figure(2);
    plot(resValues(1:count));
    xlabel('Iterations')
    ylabel('Residual Norm')
    title('Residual Norm across Iterations')
end

end

