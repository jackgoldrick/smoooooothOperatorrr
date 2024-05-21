function res = genComparison(cMatrix, p, err_a)
    if p == 1
        res = p1Norm(cMatrix);
        return

    end
    cMatrixPrime = cMatrix';
    hMatrix = cMatrixPrime * cMatrix;

    
    loc = colMaxP(cMatrix, p);
    vNow = hMatrix(:, loc);

    clear hMatrix

    q = 1 / (1 - 1/p);
    error = 1;
    res = 0;

    if p < 2  && p > 1
        temp = p;
        p = q;
        q = temp;

        temp = cMatrix;
        cMatrix = cMatrixPrime;
        cMatrixPrime = temp;
        if p == 2
            cMatrixPrime = transpose(cMatrixPrime);
            cMatrix = transpose(cMatrix);
        end

    end 

    oldGuess = vectorPNorm(vNow, q);
    vNow = vNow ./ oldGuess;
    % if p == 2
    %         cMatrixPrime = transpose(cMatrixPrime);
    %         cMatrix = transpose(cMatrix);
    % end
    
    while (error > err_a)
        
            
        vNext = cMatrixPrime * vNow;

        guess = vectorPNorm(vNext, q);

        vNextDualNormed =  vNext ./ guess;

        z = cMatrix * vNextDualNormed;

        error = abs((guess - oldGuess) / guess);

        % if Cnorm.vectorPNorm(z,q) <= z' * vNow| | error < err_a
        if error < err_a
            break;

        end

        vNow =  z ./ vectorPNorm(z, p);
        oldGuess = guess;


        
    end 
    
    res = guess;
    


end