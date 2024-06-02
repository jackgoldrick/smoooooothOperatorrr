function [res, vMax] = genComparisonSlow(cMatrix, p, err_a, sMax)

% rowSize = length(cMatrix(:,1));



if p == 1
    [res, m] = p1(cMatrix);
    vMax = cMatrix(m, :);
    return
end
cMatrixPrime = cMatrix';


q = 1 / (1 - 1/p);

error = 1; %#ok<NASGU>
res = 0; %#ok<NASGU>

if p > 2 
    temp = p;
    p = q;
    q = temp;

    temp = cMatrix;
    cMatrix = cMatrixPrime;
    cMatrixPrime = temp;


end

colSize = length(cMatrix(1,:));
x = randn(colSize, sMax,"like",1i);
vMax = ones(colSize, 1);
xNorm = vectorPNorm(x, p);
x = x ./ xNorm;
oldGuess = 0;

while (true)
    
    y = cMatrix * x;
    
    yDual = (abs(y).^ (p-1)) .* (sign(y));
    
    yDual =  yDual ./ vectorPNorm(yDual, q);
    

    [guess, index] = max(vectorPNorm(y, p));
    
    if abs(guess - oldGuess) < err_a
        vMax = y(:, index);
        break;
        
    end

    z = cMatrixPrime * yDual;
    zDual = (abs(z).^ (q-1)) .* (sign(z));
    zDual =  zDual ./ vectorPNorm(zDual, p);

    x = zDual;
    oldGuess = guess;

    if isnan(guess)
        fprintf("Floating Point Error \n")
        res = 0;
        return
    end 
end

res = guess;

end 
