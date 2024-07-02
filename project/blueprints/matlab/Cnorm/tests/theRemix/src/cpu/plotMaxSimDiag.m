function proportions = plotMaxSimDiag(N, dp, pMax, sMax, invertableMatrix, DSBMaxesDesired, bound)

    programTime = tic;

    if isempty(N) %#ok<BDLOG>
        cMatrix = [2 1; 1 2];
        N = 2;
    elseif (length(N) - 1) %#ok<BDLOG>
        cMatrix = complex(randn(N(1), N(2)), randn(N(1), N(2)));
    else 
         cMatrix = complex(randn(N, N), randn(N, N));
    end
           

    
    p = 1:dp:pMax;
    sizeP = length(p);
 
    norms = zeros(1, sizeP);
    norms2 = zeros(1, sizeP);

    %% Computes the P->P and/or P->R depending on pq flag
     % 'g' - computes P->R with O(sizeP * (sizeR + 1)) total operator norms calculated with graph
     % 'q' - computes P->R for ONE value of R ONLY
     % 'g' - computes P->P and P->R with O(sizeP * sizeR) total operator norms calculated without graph
     % 'p' - computes P->P ONLY 

    correctNorms = zeros(1, sizeP);
    proportions = zeros(1, sizeP-1, 2);
    
    if invertableMatrix == 0
        if (length(N)-1)
            N = N(2);
        end
        invertableMatrix = hadamard(N) ./ sqrt(N);
    end   

    %Needs to calculate once for first DSBMaxes set
    [norms(1), correctNorms(1), DSBMaxes] = maxSimDiag(cMatrix, invertableMatrix , p(1), 1e-1, sMax, 1e-7, DSBMaxesDesired);
    [norms2(1)] = maxSimDiag(cMatrix, invertableMatrix, p(1), 1e-1, sMax + DSBMaxesDesired, 1e-7);
    
    for j = 2:sizeP
        [norms(j), correctNorms(j), DSBMaxes, proportions(1, j-1, 1)] = maxSimDiag(cMatrix, invertableMatrix , p(j), 1e-1, sMax + min(DSBMaxesDesired, 0), 1e-7, abs(DSBMaxesDesired), DSBMaxes);
        [norms2(j), ~, proportions(1, j-1, 2)] = maxSimDiag(cMatrix, invertableMatrix, p(j), 1e-1, sMax + DSBMaxesDesired, 1e-7);
        if p(j) == 2.5
            fprintf('\nStop Here\n');
        end    
    end
    
    %% P->P Plotting from gp flag
    figure 1;

    hold on;
    plot(p, norms, '-b', 'LineWidth', 1);
    plot(p, norms2, '-c', 'LineWidth', 1);
    plot(p, correctNorms, '-r');
    ylabel('Norms');

    xlabel('p');

    title(('MaxSimDiag with BMax guesses - usefulness test'));
    
    legend("Gradient Algorithm", "pPower Algorithm");
    hold off

    figure 2;

    histogram(proportions(:,:,1), [0, .8, 1.2, 2, 4, 8]);
    hold off;

    figure 3;

    histogram(proportions(:,:,2), [0, .8, 1.2, 2, 4, 8]);
    hold off;

    programTime = toc;

    fprintf('\nTook %d \n', programTime);
   
end 