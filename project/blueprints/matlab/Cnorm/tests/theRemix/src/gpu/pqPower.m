%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%
    %%                                                                                                                         %%
        %% Author: @jackgoldrick                                                                                           %%
        %% Repository: https://github.com/jackgoldrick/smoooooothOperatorrr                                                %%
        %%                                                                                                                 %%
        %%                                                                                                                 %%
        %% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%                                                                                            
        %% 'cMatrix' - MxN Operator that the power method will calculate the norm for                                      %%
        %% 'p' - p value for the norm, this is p commonly seen in the norm subscript p->q                                  %%
        %% 'r' - r value for the norm, this is q commonly seen in the subscript p->q                                       %%
        %% 'err_a' - minimum iterative convergence tolerence                                                               %%
        %% 'sMax' - If random guesses are needed to start the algorithm, this determines the number of parallel guesses    %%
        %% 'vMax4' - a single user defined guess that can be concat to the parallel guess matrix when provided             %%
    %%                                                                                                                         %%
%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%

function [res, vMax] = pqPower(cMatrix, p, r, err_a, sMax, vMax4)

% rowSize = length(cMatrix(:,1));



if p == 1
    res = max(vectorPNorm(cMatrix, r)); 
    vMax = [];
    return
end

colSize = length(cMatrix(1,:));
vMax = zeros(colSize, 1,'gpuArray');
cMatrixPrime = cMatrix';

q = 1 / (1 - 1/p);
s = 1 / (1 - 1 ./ r);


error = 1; %#ok<NASGU>
res = 0; %#ok<NASGU>

if nargin == 6
    if isempty(vMax4)
        x = randn(colSize, sMax,"like",1i,'gpuArray');
    else
        f = randn(colSize, sMax,"like",1i,'gpuArray');
        x = [f,vMax4];
        clear f
    end

else 

    x = randn(colSize, sMax,"like",1i,'gpuArray');

end

xNorm = vectorPNorm(x, p);

x = x ./ xNorm;
oldGuess = 0;

while (true)
    
    y = cMatrix * x;
        
    yDual = (abs(y).^ (r-1)) .* (sign(y));
    
    yDual =  yDual ./ vectorPNorm(yDual,s);
    
    z = cMatrixPrime * yDual;
    
    [guess, index] = max(vectorPNorm(y, r));
    
    if abs(guess - oldGuess) < err_a
        vMax = x(:, index);
        break;
        
    end

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
