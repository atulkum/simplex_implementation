function E = createIncidenceGraph(fileName)

    fileID = fopen(fileName);
    hdr = fscanf(fileID,'p %d %d\n')

    E=sparse(zeros(hdr(2), hdr(1)));

    eI = fscanf(fileID,'%d %d\n',[2,hdr(2)] );
        
    for i = 1:hdr(2)
       E(i, eI(1,i)) = 1;
       E(i, eI(2,i)) = 1;
    end
    
    fclose(fileID);
end