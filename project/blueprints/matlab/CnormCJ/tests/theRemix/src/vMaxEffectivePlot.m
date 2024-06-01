
function vMaxEffectivePlot(N, pMax, dp, graphs, sMax, type)
    %% Test Matricies (type)
      % 'a' - 2x4 of 1s on the Diagnol 
      % 'b' - 
      % 'c' - MxN Random (Normal Distribution) Complex Matrix
      % 'd' - 3x3 Alonso's example that is elementwise divided by -3
      % 'e' - 2x2 1234 example from the slides
      % 'f' - 
      % 'g' - 4x4 example with real primes as entries
      % 'h' - NxN Hadamard Matrix
      % 'i' - NxN Identity Matrix 
      % 'j' - 4x4 Symmetric Matrix with zeros and Primes
      % 'k' - 
      % 'l' - 
      % 'm' - 3x3 Alonso's example that is elementwise divided by sqrt(3)
      % 'n' - 
      % 'o' - 
      % 'p' - 4x2 trasnpose of case 'a'
      % 'q' - 
      % 'r' - MxN Random (Normal Distribution) Real Matrix
      % 's' - 4x2 Wilson's example with 3, 1, -1 on diagnols and off diagnols
      % 't' - 
      % 'u' - 
      % 'v' -

    switch type
        case 'a'
            cMatrix = [ 1 0 1 0; 
                        0 1 0 1];

        case 'c'
            if (length(N) - 1) %#ok<BDLOG>
                cMatrix = complex(randn(N(1), N(2)), randn(N(1),N(2)));
            else
                cMatrix = complex(randn(N, N), randn(N,N));
            end

        case 'd'
            cMatrix = [ -2 1 1; 1 -2 1; 1 1 -2];
            cMatrix = cMatrix .* (1 / -3);

        case 'e'
            cMatrix = [ 1 3; 
                        2 4];

        case 'g'
            cMatrix = [ 37 0  91  0; 
                        19 73 11 47; 
                        83 97 17  0;
                        13 32 23 5];

        case 'h'
            cMatrix = hadamard(N);
            norm(cMatrix, 1)
            
        case 'i'
            cMatrix = eye(N);
            
        case 'j'
            cMatrix = [ 37 0 83  0; 
                        0 71 11 47; 
                        83 11 17 0;
                        0  47  0 5];

        case 'm'
            cMatrix = [ 0 1 -1; -1 0 1; 1 -1 0];
            cMatrix = cMatrix .* (1 / sqrt(3));

        case 'p'
            cMatrix = transpose([ 1 0 1 0; 0 1 0 1]);

        case 'r'
            if (length(N) - 1) %#ok<BDLOG>
                cMatrix = rand(N(1), N(2));

            else
                cMatrix = randn(N, N);
            end
        
        case 's'
            cMatrix = [ 3  1; 
                        1  3; 
                        3 -1;
                        -1 3];

    end

    programTime = tic;
    p = 1:dp:pMax;
    sizeP = length(p);

    
    norms = zeros(sizeP, graphs);
    

    if type == 'a' || type == 'p'
        correctNormsMax = zeros(1, sizeP);
        correctNormsQ = zeros(1, sizeP);
        correctNormsP = zeros(1, sizeP);

    else 
        correctNorms = zeros(1, sizeP);
    end


    %% Computes the P->P and/or P->R depending on pq flag
     % 'b' - computes P->R with O(sizeP * sizeR) total operator norms calculated with graph
     % 'q' - computes P->R for ONE value of R ONLY
     % 'g' - computes P->P and P->R with O(sizeP * (sizeR + 1)) total operator norms calculated without graph
     % 'p' - computes P->P ONLY
    compTime = tic;
    fprintf('Computing... \n');

    if pq == 'b' 
        fprintf('Note: Algorithm is ran %10.1f times \n', sizeP * (sizeR + 1));

    elseif pq == 'g'
        fprintf('Note: Algorithm is ran %10.1f times \n', sizeP * sizeR);
    else 
        fprintf('Note: Algorithm is ran %10.1f times \n', sizeP);

    end 

    for j = 1:sizeP    
        if pq == 'b'
            for k =1:sizeR
                [norms3D(j, k), ~] = genComparisonPQ(cMatrix, p(j), r(k), .000000001, sMax);
            end
            [norms(j), ~] = genComparison(cMatrix, p(j), .000000001, sMax);
        elseif pq == 'g'
            for k =1:sizeR
                [norms3D(j, k), ~] = genComparisonPQ(cMatrix, p(j), r(k), .000000001, sMax);
            end
        elseif pq == 'q'
            [norms(j), ~] = genComparisonPQ(cMatrix, p(j), r, .000000001, sMax);
        else      
            if ipp == 0
               [norms(j), ~] = genComparison(cMatrix, p(j), .000000001, sMax, 0);
           elseif ipp == 1
               if j == 1
                   [norms(j), vMax(:,j)] = genComparison(cMatrix, p(j), .000000001, sMax, 0);
               else
                   [norms(j), vMax(:,j)] = genComparison(cMatrix, p(j), .000000001, sMax, vMax(:,j-1));
               end
           else
               fprintf('trying %d \n', j);
               [norms(j), ~] = genComparison(cMatrix, p(j), .000000001, sMax, 0);
               fprintf('1\n');
               if j== 1
                   [norms2(j), vMax(:,j)] = genComparison(cMatrix, p(j), .000000001, 1, 0);
               else
                   [norms2(j), vMax(:,j)] = genComparison(cMatrix, p(j), .000000001, 1, vMax(:,j-1));
               end
           end
        end
    end
    fprintf('Done! \n');
    toc(compTime)
    
    
        
    if gp == 'y'
        fprintf('\nCollecting Exact Values... ')
        if type == 'h'
            for j = 1:sizeP
                q = 1 / (1 - 1 / p(j));
                correctNorms(j) = max(N ^ (1/p(j)), N ^ (1 / q));
            end 

        elseif type == 'a' || type == 'p'
            for j = 1:sizeP
                q = 1 / (1 - 1 / p(j));
                correctNormsMax(j) = max(2 ^ (1 / q), 2 ^ (1 / p(j)));
                correctNormsQ(j) = 2 ^ (1 / q);
                correctNormsP(j) = 2 ^ (1 / p(j));
            end

        elseif p(j) <= 2 && type == 'm'
            for j = 1:sizeP
                correctNorms(j) = ((1 + 2 ^(p(j) - 1)) ^ (1/(p(j)))) / sqrt(3);
            end

        elseif p(j) > 2 && type == 'm'
            for j = 1:sizeP
                q = 1 / (1 - 1 / p(j));
                correctNorms(j) = ((1 + 2 ^(q - 1)) ^ (1/(q))) / sqrt(3);
            end
        end 
        
        fprintf("Done! \n");
    end
    
    %% P->P Plotting from gp flag
    if gp == 'y'
        plotTime = tic;
        fprintf('\nPlotting P->P...')
        figure
            hold on
                xlabel("p");
                ylabel("Operator Norm Value");
                title("Operator Norm vs p");
                
                if type == 'a' || type == 'p'
                    plot(p, norms, '--b', 'LineWidth', 2);
                    plot(p, correctNormsQ, ':r', 'LineWidth', 1);
                    plot(p, correctNormsP, '-.', 'Color', "#D95319", "LineWidth", 1); 
                    plot(p, correctNormsMax, '-g', 'LineWidth', 1);
                    legend("Our Method", "Exact Q", "Exact P", "Exact Max");
                elseif type == 'i'
                    plot(p, norms, '--b', 'LineWidth', 1);
                    plot(p, N .^ (1/r  - 1./p), '-r');
                    legend("Our Method", "Exact Value");
                elseif type == 'd'
                    plot(p, norms, '--b', 'LineWidth', 1);
                    plot(p, 1, '-r');
                    legend("Our Method", "Lower Bound");
                else 
                    plot(p, norms, '--b', 'LineWidth', 1);
                    plot(p, correctNorms, '-r');
                    if ipp == 2
                        plot(p, norms2, '-y', 'LineWidth', 1);
                    end
                    legend("Our Method", "Exact Value");
                end
            hold off
        fprintf("Done! \n");
        toc(plotTime)
    end
    

    %% 3D Plot 
    if pq == 'g' 
        surfPlotTime = tic;
        fprintf('\nPlotting P->R Surface... \n');
        
        figure
            hold on
                xlabel("r");
                ylabel("p");
                zlabel("P->R Operator Norm Value");
                title("Operator Norm vs r,p");
                surf(r, p, norms3D);
                legend("Our Method");
            hold off
            
        fprintf("Done! \n ");
        toc(surfPlotTime)
    end
    
    fprintf('\nProgram Runtime:');
    toc(programTime)

end 