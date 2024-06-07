%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%
    %%                                                                                                                         %%
        %% Author: @jackgoldrick @wiwu2390                                                                                 %%
        %% Repository: https://github.com/jackgoldrick/smoooooothOperatorrr                                                %%
        %%                                                                                                                 %%
        %%                                                                                                                 %%
        %% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%                                                                                            
        %% 'N' - size of NxN matrix. If N is a vector, then for any MxN matrix, N := [M N];                                %%
        %% 'dp' - step size of the value of P                                                                              %%
        %% 'pMax' - max value of P to compute                                                                              %%
        %% 'sMaxDiag' - The number of random guess vectors used to compute each P->P or P->R                                   %%
        %% 'type' - The type of matrix to compute                                                                          %%
        %% 'gp' - Graphs the P->P depending on this flag                                                                   %%
        %%     'y' - yes P->P Graph                                                                                        %%
        %%     'n' - no P->P Graph                                                                                         %%
        %% 'pq' - Computes the P->P and/or P->R depending on this flag                                                     %%
        %%     'g' - computes P->R with O(sizeP * (sizeR + 1)) total operator norms calculated with graph                  %%
        %%     'q' - computes P->R for ONE value of q ONLY                                                                 %%
        %%     'b' - computes P->P and P->Q with O(sizeP * sizeR) total operator norms calculated without graph            %%
        %%     'p' - computes P->P ONLY                                                                                    %%
        %% 'rv' - values of R to compute. If vector => rv := [rMin rMax]. If scalar => alg runs against one value of R     %%
        %% 'dr' - step size of the value of R                                                                              %%
    %%                                                                                                                         %%
%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%

function [result, pNormStackA] = maxSimDiag(diagonalStackA, matrixU, p, step, sMaxDiag, err_a)
    dim1 = length(matrixU(1,:));
    dim2 = length(diagonalStackA(:,1));
    errorScale = 100;
    result = 0;
    q = 1 / (1 - 1/p);
    sMaxPower = 25;
    stackA = zeros(dim1 *dim2, dim1);
    stackB = zeros(dim1, dim1 * dim2, sMaxDiag);
    matrixUinv = inv(matrixU);

    for i=0:dim2 - 1
    
        stackA(((i  * dim1) + 1):((i+1) * dim1), :) = matrixU *  diag(diagonalStackA(i + 1, :)) * matrixUinv; %#ok<*MINV>
    end
    
    
    
    [pNormStackA, ~] = pPower(stackA, p, err_a, sMaxPower);


    v = zeros(dim1, 1, sMaxDiag); %#ok<PREALL>


        
    
    diagonalStackB = complex(randn(dim2, 1, dim1, sMaxDiag), randn(dim2, 1, dim1, sMaxDiag));
    
    
    
    w = complex(randn(dim1, 1, sMaxDiag), randn(dim1, 1,sMaxDiag));
    
    
    
    pNormBA = 0; %#ok<NASGU>
    
    oldGuess = 0;
    vMax = zeros(dim1 * dim2, sMaxDiag);
    
    
    while (true)
        
    
        stackB = pagemtimes(matrixU,  diagonalStackB .* matrixUinv);
        stackB = reshape(stackB, dim2 * dim1, dim1, sMaxDiag);
         
        [normStackB, vMax] = pPower(stackB, p, err_a, sMaxPower, vMax);
        diagonalStackB = diagonalStackB ./ normStackB;
        diagSum = sum((diagonalStackA .* diagonalStackB), 1);
    
        vDual = ((w' * matrixU) .* diagSum) * matrixUinv;
        v = dual(vDual.', q);
        % vNorm = vectorPNorm(v.', p);
        wDual =  matrixU * (diagSum.' .* (matrixUinv * v));
        w = dual(wDual, p);
        % wNorm = vectorPNorm(w, q);
        wPrime_U = w' * matrixU;
        matrixUinv_v = matrixUinv * v;
    
        guess = wPrime_U * (diagSum.' .* matrixUinv_v);
        guess = abs(guess);
        gradient = zeros(dim2, dim1);
    
        for i=1:dim2
    
            gradient(i, :) = (wPrime_U) .* diagonalStackA(i,:) .* (matrixUinv_v).';
        end
    
        if (guess > pNormStackA + errorScale * err_a)
            fprintf('Matrix Nrom inequality violated \n');
            fprintf('  diff = %d \n', pNormStackA - guess);
            break;
        end
    
        if abs(guess - pNormStackA) <=err_a || (guess - oldGuess) <= err_a
            break;
        end
    
        diagonalStackB = diagonalStackB +  step .* gradient;
    
        oldGuess = guess;
    
    
    
    end

    result = max(result, guess);


    % pNormStackA;
    diff = pNormStackA - result; %#ok<NOPRT>

    if diff < - errorScale * err_a
        fprintf('\n  Something went all Fucky-Whucky \n')
    end
        

end