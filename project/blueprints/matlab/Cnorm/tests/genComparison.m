function res = genComparison(cMatrix, p, err_a)
    rowSize = length(cMatrix(:,1));
    colSize = length(cMatrix(1,:));
    % I currently have zero clue why this fixes things
    
    if p == 1 
        if colSize >= rowSize 
            if isreal(cMatrix(1,1))
                res = pInf(cMatrix)
                return
            else 
                cMatrixPrime = cMatrix';
                conjTransMatrix = cMatrixPrime * cMatrix;

                res = pInf(conjTransMatrix)
                % res = p1(conjTransMatrix);
                % res = pInf(cMatrix);
                % res = p1(cMatrix);
                clear conjTransMatrix
                return
            end

        else
          
            res = p1(cMatrix);
            return
             
        end
        % res = p1(cMatrix);
        % return
    end
    cMatrixPrime = cMatrix';


    


    q = 1 / (1 - 1/p);
    error = 1;
    
    % commenting this block made the 2 ^ 1/q function exact with A and W's example
    % If removed Alonso's Example breaks
    %  This handles  the p-q dual between 1-2
    if p < 2  && p > 1
        temp = p;
        p = q;
        q = temp;

        temp = cMatrix;
        cMatrix = cMatrixPrime;
        cMatrixPrime = temp;


    end 


    
    loc = colMaxP(cMatrixPrime, p);
    vNow = cMatrix(:, loc);

    oldGuess = vectorPNorm(vNow, q);
    vNow = vNow ./ oldGuess;
    
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