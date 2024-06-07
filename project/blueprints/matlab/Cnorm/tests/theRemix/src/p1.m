function [res, ind] = p1(matrix)

    modMatrix = abs(matrix);
    normColumn =  sum(modMatrix, 1);

    [res, ind] = max(normColumn,[], 2);

  

end 