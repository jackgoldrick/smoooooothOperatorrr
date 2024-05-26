function res = p1(matrix)
    oneVector = ones(1,length(matrix(:,1)));
    % modMatrix = sqrt((real(cMatrix)) .^ 2 + (imag(cMatrix)) .^ 2);
    % on avg saves a bit of time compared to above line
    modMatrix = abs(matrix);
    normColumn =  oneVector * modMatrix;

    [res, ~] = max(normColumn);

end 