function [res, ind] = pInf(matrix)
    modMatrix = abs(matrix);
    normColumn =  sum(modMatrix, 2);

    [res, ind] = max(normColumn,[], 1);
    
end