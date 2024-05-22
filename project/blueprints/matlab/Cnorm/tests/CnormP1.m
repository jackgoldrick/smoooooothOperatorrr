function res = CnormP1(N)
    tic
    % creates random NxN C^d matrix
    cMatrix = complex(rand(N, N), rand(N,N));


    % creates 
    oneVector = ones(1,N);
    % modMatrix = sqrt((real(cMatrix)) .^ 2 + (imag(cMatrix)) .^ 2);
    % on avg saves a bit of time compared to above line
    modMatrix = abs(cMatrix);
    normColumn =  oneVector * modMatrix;

    [res, ~] = max(normColumn);
    
    toc

    norm(cMatrix, 1)
end