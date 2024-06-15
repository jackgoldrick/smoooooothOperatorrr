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

function [result, pNormStackA, DSBMaxes] = maxSimDiag(diagonalStackA, matrixU, p, step, sMaxDiag, err_a, DSBMaxDesired, diagonalStackBGuess)
DSAOriginal = diagonalStackA;
dims = size(diagonalStackA);
dim1 = dims(2);
dim2 = dims(1);
matrixUinv = inv(matrixU);
diagonalStackA = reshape(diagonalStackA.', dim1, 1, dim2);

% reshapedU = reshape(matrixU, dim1, dim1, 1);
% reshapedUinv = reshape(matrixUinv, dim1, dim1, 1);

errorScale = 100;
result = 0;
q = 1 / (1 - 1/p);
sMaxPower = 25;

diagonalStackB = complex(randn(dim1, 1, dim2, sMaxDiag), randn(dim1, 1, dim2, sMaxDiag));

if nargin == 8
    diagonalStackB = cat(4, diagonalStackB, diagonalStackBGuess);
    dimsTemp = size(diagonalStackB);
    sMaxDiag = dimsTemp(4);
    clear dimsTemp;
end

vMax = zeros(dim1 * dim2, 1, sMaxDiag);
step = step * ones(1,1,1,sMaxDiag);

stackA = zeros(dim1 *dim2, dim1);
stackB = zeros(dim1, dim1 * dim2, sMaxDiag);


stackA = pagemtimes(matrixU, diagonalStackA .* matrixUinv);
stackA = reshape(permute(stackA, [1 3 2]), dim2 * dim1, dim1, 1);



[pNormStackA, ~] = pPower(stackA, p, err_a, sMaxPower);


v = zeros(dim1, 1, 1, sMaxDiag); %#ok<PREALL>


w = complex(randn(dim1, 1, 1, sMaxDiag), randn(dim1, 1, 1, sMaxDiag));


pNormBA = 0; %#ok<NASGU>

oldGuess = 0;

converged = zeros(1, 1, 1, sMaxDiag);

while (true)

    % reshape(matrixUinv,1, dim1, dim1,1)
    stackB = pagemtimes(matrixU, diagonalStackB .* matrixUinv);
    stackB = reshape(stackB, dim1, dim2 * dim1, sMaxDiag);
    
    flatConverged = reshape(converged, 1, []).';
    normStackB = ones(1, 1, sMaxDiag);
    [normStackB(:,:,~flatConverged), vMax(:,:,~flatConverged)] = pPower(stackB(:,:,~flatConverged), p, err_a, sMaxPower, vMax(:,:,~flatConverged));
    diagonalStackB = diagonalStackB ./ reshape(normStackB, 1, 1, 1, sMaxDiag);
    % diagSum = reshape(sum((diagonalStackA .* diagonalStackB), 3), dim1, 1, sMaxDiag);
    diagSum = sum((diagonalStackA .* diagonalStackB), 3);
    vDual = pagemtimes((pagemtimes(pagectranspose(w), matrixU) .* pagetranspose(diagSum)), matrixUinv);
    v = dual(pagectranspose(vDual), q);
    % vNorm = vectorPNorm(v.', p);
    wDual =  pagemtimes(matrixU, diagSum .* pagemtimes(matrixUinv, v));
    w = dual(wDual, p);
    % wNorm = vectorPNorm(w, q);
    wPrime_U = pagemtimes(pagectranspose(w), matrixU);
    matrixUinv_v = pagemtimes(matrixUinv, v);

    guess = pagemtimes(wPrime_U, (diagSum .* matrixUinv_v));
    guess = abs(guess);

    [maxGuess, ~] = max(guess, [], "all");
    if (maxGuess > pNormStackA + errorScale * err_a)
        fprintf('Matrix Nrom inequality violated \n');
        pNormStackA
        guess
        %fprintf('  diff = %d \n', pNormStackA - guess);
        break;
    end

    converged = converged | abs(guess - pNormStackA) <= err_a | guess - oldGuess <= err_a;

    %if sum(converged) > 0.7 * sMaxDiag
    if all(converged)
        break;
    end

    gradient = pagetranspose(wPrime_U) .* diagonalStackA .* (matrixUinv_v);
    step(converged) = 0;
    diagonalStackB = diagonalStackB +  step .* gradient;

    oldGuess = guess;



end


result = maxGuess;
[~, indexMaxMultiple] = maxk(guess, DSBMaxDesired);
DSBMaxes = diagonalStackB(:, :, :, indexMaxMultiple);


% pNormStackA;
diff = pNormStackA - result; %#ok<NOPRT>

if min(diff,[],"all") < - errorScale * err_a
    fprintf('\n  Something went all Fucky-Whucky \n')
end


end
