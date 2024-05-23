function res = genComparison(cMatrix, p, err_a)
rowSize = length(cMatrix(:,1));
colSize = length(cMatrix(1,:));
% I currently have zero clue why this fixes things

if p == 1
    % if colSize >= rowSize
    %     if isreal(cMatrix(1,1))
    %         res = pInf(cMatrix) %#ok<NOPRT>
    %         lab = norm(cMatrix, inf) %#ok<NOPRT,NASGU>
    %         % res = p1(cMatrix)
    %         % lab = norm(cMatrix, 1)
    %         condition = cond(cMatrix) %#ok<NOPRT,NASGU>
    %         return
    %     else
    %         % cMatrixPrime = cMatrix';
    %         % conjTransMatrix = cMatrixPrime * cMatrix;
            
    %         % res = pInf(conjTransMatrix) %#ok<NOPRT>
    %         % res = p1(conjTransMatrix);
    %         % res = pInf(cMatrix);
    %         res = p1(cMatrix) %#ok<NOPRT>
    %         clear conjTransMatrix
    %         return
    %     end
        
    % else
        
    %     res = p1(cMatrix);
    %     return
        
    % end
    res = p1(cMatrix);
    return
end
cMatrixPrime = cMatrix';





q = 1 / (1 - 1/p);
error = 1;

% commenting this block made the 2 ^ 1/q function exact with A and W's example
% If removed Alonso's Example breaks
%  This handles  the p-q dual between 1-2
% if p < 2  && p > 1
%     temp = p;
%     p = q;
%     q = temp;
% 
%     temp = cMatrix;
%     cMatrix = cMatrixPrime;
%     cMatrixPrime = temp;
% 
% 
% end



%loc = colMaxP(cMatrixPrime, p);
x = cMatrix(1, :)';

xNorm = vectorPNorm(x, p);
x = x ./ xNorm;
oldGuess = 0;

%while (error > err_a)
while (true)
    
    y = cMatrix * x;
    
    % vNext = vNext ./ vectorPNorm(vNext, inf);
    
    yDual = (abs(y).^ (p-1)) .* (sign(y));
    
    yDual =  yDual ./ vectorPNorm(yDual, q);
    
    z = cMatrixPrime * yDual;
    
    %error = abs((guess - oldGuess) / guess);
    guess = vectorPNorm(y, p);
    
    %if vectorPNorm(z,q) <= z' * x + 1e-6%| | error < err_a
    if abs(guess - oldGuess) < err_a
        break;
        
    end
    % These lines can be commented out
    %z = z ./ vectorPNorm(z, inf);
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