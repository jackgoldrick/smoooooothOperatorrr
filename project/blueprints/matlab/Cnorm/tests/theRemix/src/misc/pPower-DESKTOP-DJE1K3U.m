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

dims = size(cMatrix);
colSize = dims(2);
if length(dims) == 3
    pageDim = dims(3);

else
    pageDim = 1;

end 

if p == 1
    [res, index] = p1(cMatrix);
    vMax = zeros(colSize, 1, pageDim);
    vMax(index == 1:colSize) = 1;
    return
end

cMatrixPrime = pagectranspose(cMatrix);


q = 1 / (1 - 1/p);

error = 1; %#ok<NASGU>
res = 0; %#ok<NASGU>







if nargin == 5
    if isempty(vMax3)
        x = randn(colSize, sMax, pageDim, "like",1i);
    else
        f = randn(colSize, sMax, pageDim, "like",1i);
        x = [f,vMax3];
        sMax= sMax + 1;
        clear f
    end

else 

    x = randn(colSize, sMax, pageDim, "like",1i);

end



xNorm = vectorPNorm(x, p);
x = x ./ xNorm;
oldGuess = 0;

while (true)
    
    y = pagemtimes(cMatrix, x);
    
    yDual = dual(y, p);

    [guess, index] = max(vectorPNorm(y, p), [], 2);    
    
    if max(abs(guess - oldGuess)) < err_a

        vMax(:, 1, :) = x(:, index == 1:sMax);
        break;
    end

    z = pagemtimes(cMatrixPrime, yDual);
    
    x = dual(z, q);
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
