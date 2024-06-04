%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%
        %%                                                                                                                 %%
        %%                                                                                                                 %%
        %% Author: @jackgoldrick                                                                                           %%
        %% Repository: https://github.com/jackgoldrick/smoooooothOperatorrr                                                %%
        %%                                                                                                                 %%
        %%                                                                                                                 %%
        %% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%                                                                                           
        %%                                                                                                                 %%
        %% 'cMatrix' - MxN Operator that the power method will calculate the norm for                                      %%
        %% 'p' - p value for the norm, this is p commonly seen in the norm subscript p->q                                  %%
        %% 'err_a' - minimum iterative convergence tolerence                                                               %%
        %% 'sMax' - If random guesses are needed to start the algorithm, this determines the number of parallel guesses    %%
        %% 'vMax3' - a single user defined guess that can be concat to the parallel guess matrix when provided             %%
        %%                                                                                                                 %%
%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%
function [res, vMax] = pPower(cMatrix, p, err_a, sMax, vMax3)

colSize = length(cMatrix(1,:));

if p == 1
    [res, index] = p1(cMatrix);
    vMax = zeros(colSize, 1);
    vMax(index, 1) = 1;
    return
end

cMatrixPrime = cMatrix';


q = 1 / (1 - 1/p);

error = 1; %#ok<NASGU>
res = 0; %#ok<NASGU>







if nargin == 5
    if isempty(vMax3)
        x = randn(colSize, sMax,"like",1i);
    else
        f = randn(colSize, sMax,"like",1i);
        x = [f,vMax3];
        clear f
    end

else 

    x = randn(colSize, sMax,"like",1i);

end



xNorm = vectorPNorm(x, p);
x = x ./ xNorm;
oldGuess = 0;

while (true)
    
    y = cMatrix * x;
    
    yDual = (abs(y).^ (p-1)) .* (sign(y));
    
    yDual =  yDual ./ vectorPNorm(yDual, q);
    

   [guess, index] = max(vectorPNorm(y, p));
    
    
    if abs(guess - oldGuess) < err_a
        vMax = x(:, index);
        break;
    end

    z = cMatrixPrime * yDual;
    zDual = (abs(z).^ (q-1)) .* (sign(z));
    zDual =  zDual ./ vectorPNorm(zDual, p);

    x = zDual;
    oldGuess = guess;

    if isnan(guess)
        fprintf("Floating Point Error \n")
        vMax = [];
        res = 0;
        return
    end 
end

res = guess;

end 
