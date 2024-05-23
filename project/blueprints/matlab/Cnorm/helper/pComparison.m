function [m, o, diff, N, p]= pComparison(N, p)
    addpath ../src/
        
    cMatrix = Cnorm(N);
    if p ==2
        tic
        m = norm(cMatrix.cMatrix, p)
        toc
    end

    o = real(Cnorm.pPower(cMatrix, p, .0000001))
    if p == 2 
        diff = abs(m - o)
    end 



end