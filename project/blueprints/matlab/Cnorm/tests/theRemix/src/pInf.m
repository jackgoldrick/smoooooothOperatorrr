function [res, ind] = pInf(Matrix)
    oneVector = ones(length(Matrix(1,:)), 1);
    % modMatrix = sqrt((real(cMatrix)) .^ 2 + (imag(cMatrix)) .^ 2);
    % on avg saves a bit of time compared to above line
    modMatrix = abs(Matrix);
    normColumn =  modMatrix * oneVector;
   
    [res, ind] = max(normColumn);

end