function [res, vMax] = genComparison(cMatrix, p, err_a, sMax)

    % rowSize = length(cMatrix(:,1));
    colSize = length(cMatrix(1,:));
    vMax = zeros(colSize, 1);


    if p == 1
        res = p1(cMatrix) %#ok<NOPRT>
        return
    end
    cMatrixPrime = cMatrix';


    q = 1 / (1 - 1/p);
    error = 1;


    % if p > 2 
    %     temp = p;
    %     p = q;
    %     q = temp;

    %     temp = cMatrix;
    %     cMatrix = cMatrixPrime;
    %     cMatrixPrime = temp;


    % end



    res = 0;
    for j = 1:sMax
        x = randn(colSize, 1,"like",1i);


        xNorm = vectorPNorm(x, p);
        x = x ./ xNorm;
        oldGuess = 0;

        %while (error > err_a)
        while (true)
            
            y = cMatrix * x;
            
            % y = y ./ vectorPNorm(y, inf);
            
            yDual = (abs(y).^ (p-1)) .* (sign(y));
            
            yDual =  yDual ./ vectorPNorm(yDual, q);
            
            z = cMatrixPrime * yDual;
            
            %error = abs((guess - oldGuess));
            guess = vectorPNorm(y, p);
            
            %if error < err_a
            if abs(guess - oldGuess) < err_a
                vMax = y;
                break;
                
            end
            % 
            % z = z ./ vectorPNorm(z, inf);
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
        res =  guess;
    end 
    % res = guess;



end