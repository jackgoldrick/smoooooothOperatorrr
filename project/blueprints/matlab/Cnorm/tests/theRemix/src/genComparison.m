function [res, vMax] = genComparison(cMatrix, p, err_a, sMax)

% rowSize = length(cMatrix(:,1));
colSize = length(cMatrix(1,:));
vMax = ones(colSize, 1);

if p == 1
    res = p1(cMatrix) %#ok<NOPRT>
    return
end
cMatrixPrime = cMatrix';


q = 1 / (1 - 1/p);
error = 1;





% loc = colMaxP(cMatrixPrime, p);
% x1 = cMatrix(loc, :)';

% loc = colMaxP(cMatrixPrime, q);
% x3 = cMatrixPrime(loc, :)';

res = 0;


% for j = 1:sMax
    x = randn(colSize, sMax,"like",1i);
    % x = [x, start];

    
    xNorm = vectorPNorm(x, p);
    x = x ./ xNorm;
    oldGuess = 0;

    %while (error > err_a)
    while (true)
        
        y = cMatrix * x;
        
        
        yDual = (abs(y).^ (p-1)) .* (sign(y));
        
        yDual =  yDual ./ vectorPNorm(yDual, q);
        
        z = cMatrixPrime * yDual;
        
        %error = abs((guess - oldGuess));
        guess = max(vectorPNorm(y, p));
        
        %if error < err_a
        if abs(guess - oldGuess) < err_a
            vMax = y;
            break;
            
        end

        zDual = (abs(z).^ (q-1)) .* (sign(z));
        zDual =  zDual ./ vectorPNorm(zDual, p);

        x = zDual;
        oldGuess = guess;


        if isnan(guess)
            fprintf("Fuck \n")
            res = 0;
            return
        end 
    end
    
    res = guess;

end 



