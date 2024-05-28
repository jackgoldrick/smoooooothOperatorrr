function [res, vMax] = genComparisonPQ(cMatrix, p, r, err_a, sMax)

% rowSize = length(cMatrix(:,1));
colSize = length(cMatrix(1,:));
vMax = zeros(colSize, 1);


if p == 1
    res = max(vectorPNorm(cMatrix, r)) %#ok<NOPRT>
    return
end

cMatrixPrime = cMatrix';


% if p > 2 
%     temp = p;
%     p = r;
%     r = temp;

%     temp = cMatrix;
%     cMatrix = cMatrixPrime;
%     cMatrixPrime = temp;
%     p = 1 / (1 - 1/p);
%     r = 1 / (1 - 1/r);


% end

q = 1 / (1 - 1/p);
s = 1 / (1 - 1/r);


error = 1;





% loc = colMaxP(cMatrixPrime, p);
% x1 = cMatrix(loc, :)';

% loc = colMaxP(cMatrixPrime, q);
% x3 = cMatrixPrime(loc, :)';

res = 0;


% for j = 1:sMax
    x = randn(colSize, sMax,"like",1i);
    % x = start;

    
    xNorm = vectorPNorm(x, p);
    x = x ./ xNorm;
    oldGuess = 0;

    %while (error > err_a)
    while (true)
        
        y = cMatrix * x;
        
        
        yDual = (abs(y).^ (r-1)) .* (sign(y));
        
        yDual =  yDual ./ vectorPNorm(yDual,s);
        
        z = cMatrixPrime * yDual;
        
        %error = abs((guess - oldGuess));
        guess = max(vectorPNorm(y, r));
        
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



