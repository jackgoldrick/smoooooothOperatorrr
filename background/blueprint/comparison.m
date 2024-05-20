function comparison(N)
        
    
    cMatrix = Cnorm(N);

    % tic
    % m = norm(cMatrix.cMatrix, 2)
    % toc

    o = real(Cnorm.P2power(cMatrix, .0000001))

    % diff = abs(m - o)
    

   cMatrix.cMatrix =  cMatrix.cMatrix ./ o;
   cMatrix.hMatrix  = transpose(conj(cMatrix.cMatrix)) * cMatrix.cMatrix;

   oNew = real(Cnorm.P2power(cMatrix, .0000001))

end