function [xOpt, fVal,stat] = solveVertexCover(fileName)

    E = createIncidenceGraph(fileName);
    [m n] = size(E);
    %% the problem to be encoded is 
    %% min. 1^t x s.t. E * x >= 1, x binary
    c= ones(n,1);
    vType='';
    for i = 1:n
        vType=strcat(vType,'B');
    end
    
    cType='';
    for j = 1:m
        cType=strcat(cType,'L');
    end
    
    b = ones(m,1);
    
    vType
    cType
    param.msglev =1;
    timer=tic;
    [xOpt,fVal,stat] = glpk(c,E,b,[],[],cType,vType,1,param);
    tElapsed = toc(timer);
    
    
    if (stat == 5)
       fprintf ('Solved successfully...\n'); 
       fprintf ('Optimal cover has %d vertices\n', fVal);
       for i = 1:n
        if (xOpt(i) >= 1)
               fprintf ('Vertex %d is part of the cover \n', i);
        end
       end
    end
    fprintf('Time taken: %f \n',tElapsed);
    
   
end