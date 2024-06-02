%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%
    %%                                                                                                                         %%
        %% Author: @jackgoldrick @wiwu2390                                                                                           %%
        %% Repository: https://github.com/jackgoldrick/smoooooothOperatorrr                                                %%
        %%                                                                                                                 %%
        %%                                                                                                                 %%
        %% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%                                                                                            
        %% 'N' - size of NxN matrix. If N is a vector, then for any MxN matrix, N := [M N];                                %%
        %% 'dp' - step size of the value of P                                                                              %%
        %% 'pMax' - max value of P to compute                                                                              %%
        %% 'sMax' - The number of random guess vectors used to compute each P->P or P->R                                   %%
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

function [guess, pNormStackA] = maxSimDiag(diagonalStackA, matrixU, p, step, err_a)
    dim1 = length(matrixU(1,:));
    dim2 = length(diagonalStackA(:,1));
    error = 1;
    q = 1 / (1 - 1/p);
    % if  arg(1) == 'r'
    %     if (length(N) - 1) %#ok<BDLOG>
    %         for i=1:arg(2)
    %             xMatrix = complex(randn(N(1), N(2)), randn(N(1),N(2)));
    %             yMatrix = complex(randn(N(1), N(2)), randn(N(1),N(2)));
    
    
    %     else
    %         xMatrix = complex(randn(N, N), randn(N,N));
    %         yMatrix = complex(randn(N, N), randn(N,N));
    %     end
    % end
    
    
    diagonalStackB = complex(randn(dim2, dim1), randn(dim2, dim1));
    
    v = complex(randn(dim1, 1), randn(dim1, 1));
    
    w = complex(randn(dim1, 1), randn(dim1, 1));
    
    matrixUinv = inv(matrixU);
    
    stackA = zeros(dim1 *dim2, dim1);
    stackB = zeros(dim1, dim1 * dim2);
    
    for i=0:dim2 - 1
    
        stackA(((i  * dim1) + 1):((i+1) * dim1), :) = matrixU *  diag(diagonalStackA(i + 1, :)) * matrixUinv;
    end
    
    
    
    
    [pNormStackA, vMax] = pPower(stackA, p, .000000001, 25);
    
    pNormStackA
    pNormBA = 0;
    
    oldGuess = 0;
    
    
    
    while (true)
        for i=0:dim2 - 1
    
            stackB(:, ((i  * dim1) + 1):((i+1) * dim1)) = matrixU *  diag(diagonalStackB(i + 1, :)) * matrixUinv;
        end
        [normStackB, vMax] = pPower(stackB, p, .000000001, 25, vMax);
        diagonalStackB = diagonalStackB ./ normStackB;
        diagSum = sum((diagonalStackA .* diagonalStackB), 1);
    
        vDual = ((w' * matrixU) .* diagSum) * matrixUinv;
        v = dual(vDual, q).';
        vNorm = vectorPNorm(v.', p);
        wDual =  matrixU * (diagSum.' .* (matrixUinv * v));
        w = dual(wDual, p);
        wNorm = vectorPNorm(w, q);
        wPrime_U = w' * matrixU;
        matrixUinv_v = matrixUinv * v;
    
        guess = wPrime_U * (diagSum.' .* matrixUinv_v);
        guess = abs(guess);
        gradient = zeros(dim2, dim1);
    
        for i=1:dim2
    
            gradient(i, :) = (wPrime_U) .* diagonalStackA(i,:) .* (matrixUinv_v).';
        end
    
        if (guess > pNormStackA)
            fprintf('Matrix Nrom inequality violated \n');
            break;
        end
    
        if abs(guess - pNormStackA) <=err_a || abs(oldGuess - guess) <= err_a
            break;
        end
    
        diagonalStackB = diagonalStackB +  step .* gradient;
    
        oldGuess = guess;
      
    
    
    end

    if p == 2
            diff = pNormStackA - guess

            if diff < 0
                fprintf('\n  Something went all Fucky-Whucky \n')
            end
    end

end