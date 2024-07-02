%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%
    %%                                                                                                                         %%
        %% Author: @jackgoldrick                                                                                           %%
        %% Repository: https://github.com/jackgoldrick/smoooooothOperatorrr                                                %%
        %% Tensor Networks                                                                                                 %%
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
function genCompPlot(type, N, dp, pMax, sMax, pq, gp, rv, dr, invertableMatrix, cord, bound)
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
      % 's' - Simultanious Diagnolizable Stacks
      % 't' - 
      % 'u' - user defined
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
                cMatrix = randn(N(1), N(2));

            else
                cMatrix = randn(N, N);
            end
        
        case 's'
            if isempty(N) %#ok<BDLOG>
                cMatrix = [2 1; 1 2];
                N = 2;
            elseif (length(N) - 1) %#ok<BDLOG>
                cMatrix = complex(randn(N(1), N(2)), randn(N(1), N(2)));

            else 
                 cMatrix = complex(randn(N, N), randn(N, N));
            end
           

    end

    programTime = tic;
    p = 1:dp:pMax;
    sizeP = length(p);
    if nargin >= 9
        if pq == 'g'
            r = rv(1):dr:rv(2);
            sizeR = length(r);
            norms3D = zeros(sizeP, sizeR);

        else
            norms = zeros(sizeP, 1);
            
            %For Plotting Upper Bound
            if bound == 'y'
                normsBound1 = zeros(sizeP, 1);
                normsBound2 = zeros(sizeP, 1);
                infNorm = pInf(cMatrix);
                twoNorm = norm(cMatrix, 2);
            end
            
            if ~(length(rv) - 1) %#ok<BDLOG>
                r = rv;
            else 
                r = rv(1);
            end
        end 
        sizeR = length(r);
    else
        sizeR = 1;
    end 
    

    if type == 'a' || type == 'p'
        correctNormsMax = zeros(1, sizeP);
        correctNormsQ = zeros(1, sizeP);
        correctNormsP = zeros(1, sizeP);

    else 
        correctNorms = zeros(1, sizeP);
    end


    %% Computes the P->P and/or P->R depending on pq flag
     % 'g' - computes P->R with O(sizeP * (sizeR + 1)) total operator norms calculated with graph
     % 'q' - computes P->R for ONE value of R ONLY
     % 'g' - computes P->P and P->R with O(sizeP * sizeR) total operator norms calculated without graph
     % 'p' - computes P->P ONLY 
    compTime = tic;
    fprintf('Computing... \n');

    if pq == 'g' 
        fprintf('Note: P->P and P->R Algorithms ran %10.f times \n', (sizeP * (sizeR + 1)));

    elseif pq == 'b'
        fprintf('Note: P->R Algorithm is ran %10.f times \n', sizeP * sizeR);
    else 
        fprintf('Note: P->P Algorithm is ran %10.f times \n', sizeP);

    end 
    vMax2 = zeros(length(cMatrix(1,:)), 1);
    vs = zeros (length(cMatrix(1,:)), 1, sizeP);
    for j = 1:sizeP    
        if pq == 'g'
            for k = 1:sizeR
                [norms3D(j, k), ~] = pqPower(cMatrix, p(j), r(k), .000000001, sMax);
            end
            [norms(j), vMax] = pPower(cMatrix, p(j), .000000001, sMax, vMax2);
            vMax2 = vMax;
        elseif pq == 'b'
            for k =1:sizeR
                [norms3D(j, k), ~] = pqPower(cMatrix, p(j), r(k), .000000001, sMax);
            end
        elseif pq == 'q'
            [norms(j), ~] = pqPower(cMatrix, p(j), r, .000000001, sMax);
        else
            if type == 's'

            else
                [norms(j), vMax] = pPower(cMatrix, p(j), .000000001, sMax, vMax2);
                if bound == 'y'
                    normsBound1(j) = norms(1)^(1/p(j)) * infNorm^(1-1/p(j));
                    if p(j) < 2
                        normsBound2(j) = norms(1)^(2/p(j)-1) * twoNorm^(2-2/p(j));
                    else
                        normsBound2(j) = infNorm^(1-2/p(j)) * twoNorm^(2/p(j));
                    end        
                end
                vMax2 = vMax;
                vs(:,1, j) = vMax; 
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

        elseif type == 'm'
            for j = 1:sizeP
                if p(j) <= 2 
                    correctNorms(j) = ((1 + 2 ^(p(j) - 1)) ^ (1/(p(j)))) / sqrt(3);
                else
                    q = 1 / (1 - 1 / p(j));
                    correctNorms(j) = ((1 + 2 ^(q - 1)) ^ (1/(q))) / sqrt(3);
                end
            end

        elseif type == 's'
            BMaxesDesired = floor(sMax/2); %You can comment this out as you see fit to add input arg that does this with some customization
            DSBMaxes = 0;
            if invertableMatrix ~= 0
                for j = 1:sizeP
                    [norms(j), correctNorms(j), DSBMaxes] = maxSimDiag(cMatrix, invertableMatrix , p(j), 1e-1, sMax, 1e-7);

                end 
            else 

                if (length(N) - 1) %#ok<BDLOG>

                    N = N(2);
                end

            
                DSBMaxDesired = 10;
                
                [norms(1), correctNorms(1), DSBMaxes] = maxSimDiag(cMatrix, hadamard(N) ./ sqrt(N), p(1), 1e-1, sMax, 1e-7, DSBMaxDesired);
                
                for j = 1:sizeP

                    [norms(j), correctNorms(j), DSBMaxes] = maxSimDiag(cMatrix, hadamard(N) ./ sqrt(N), p(j), 1e-1, sMax, 1e-7, DSBMaxDesired, DSBMaxes);
                    
              
                end 

            end
        end 
        
        fprintf("Done! \n");
    end
    
    %% P->P Plotting from gp flag
    if gp == 'y'
        plotTime = tic;
        fprintf('\nPlotting P->P...')
        figure1 = figure; %#ok<NASGU>
            hold on
                xlabel("p");
                ylabel("Operator Norm Value");
                title("Operator Norm vs p");
                
                if type == 'a' || type == 'p'
                    plot(p, norms, '-b', 'LineWidth', 2);
                    plot(p, correctNormsQ, ':r', 'LineWidth', 1);
                    plot(p, correctNormsP, '-.', 'Color', "#D95319", "LineWidth", 1); 
                    plot(p, correctNormsMax, '-g', 'LineWidth', 1);
                    legend("Our Method", "Exact Q", "Exact P", "Exact Max");
                elseif type == 'i'
                    plot(p, norms, '-b', 'LineWidth', 1);
                    plot(p, N .^ (1/r  - 1./p), '-r');
                    legend("Our Method", "Exact Value");
                elseif type == 'd'
                    plot(p, norms, '-b', 'LineWidth', 1);
                    plot(p, 1, '-r');
                    legend("Our Method", "Lower Bound");
                elseif type == 's' 
                    plot(p, norms, '-b', 'LineWidth', 1);
                    plot(p, correctNorms, '-r');
                    plot(p, normsBound1, '-g');
                    plot(p, normsBound2, '-c');
                    legend("Gradient Algorithm", "pPower Algorithm", "Upper Bound no 1", "Upper Bound no 2");
                else 
                    plot(p, norms, '-b', 'LineWidth', 1);
                    plot(p, correctNorms, '-r');
                    plot(p, normsBound1, '-g');
                    plot(p, normsBound2, '-c');
                    legend("Our Method", "Exact Value", "Upper Bound no 1", "Upper Bound no 2");
                end
            hold off
        fprintf("Done! \n");
        toc(plotTime)
    end
    
    % 1 -> k
    % 2 -> wn
    % 3 -> eta
  
        % model = @(x) ((abs(x(1)) .* (x(2) .^ 2)) ./ sqrt((x(2).^2 + p.^2).^2 + 4 .* (x(2) .*(x(3) .* p).^2))) - norms;

        % hold on
        %     xlabel("p");
        %     ylabel("Operator Norm Value");
        %     title("FFT of Operator Norm vs p");
        %     regNorms = lsqnonlin(model, [1 1 1]);
        %     plot(p,  model(regNorms), '-g', 'LineWidth', 1);
        % hold off

    %% 3D Plot 
    if pq == 'g' 
        surfPlotTime = tic;
        fprintf('\nPlotting P->R Surface... \n');
        
        figure2 = figure; %#ok<NASGU>
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


  if cord == 'y'
        x(1,:) = abs(vs(1, 1, :));
        y(1, :) = abs(vs(2, 1, :));
        size(x)
        figure3 = figure;
            hold on
                xlabel("x-cord");
                ylabel("y-cord");
                title("maximizing coords: abs");
                plot(x, y);
                legend("Our Method"); 
            hold off 

        xRe(1,:) = real(vs(1, 1, :));
        xIm(1,:) = imag(vs(1, 1, :));

        yRe(1,:) = real(vs(2, 1, :));
        yIm(1, :) = imag(vs(2, 1, :));
        figure4 = figure;
            hold on
                tiledlayout(2,1)

                ax1 = nexttile;
                plot(ax1, xRe, xIm);
                xlabel(ax1, "Re(x-cord)");
                ylabel(ax1, "Im(x-cord)");
                title(ax1, "maximizing x-cords: im vs re");
                

                ax2 = nexttile;
                plot(ax2, yRe, yIm);
                xlabel(ax2, "Re(y-cord)");
                ylabel(ax2, "Im(y-cord)");
                title(ax2, "maximizing y-cords: im vs re");
                
            hold off 

            figure5 = figure;
            hold on
                tiledlayout(2,1)

                ax1 = nexttile;
                plot(ax1, xRe,yRe);
                xlabel(ax1, "Re(x-cord)");
                ylabel(ax1, "Re(y-cord)");
                title(ax1, "maximizing real-cords: y v x");
                

                ax2 = nexttile;
                plot(ax2, xIm, yIm);
                xlabel(ax2, "Im(x-cord)");
                ylabel(ax2, "Im(y-cord)");
                title(ax2, "maximizing imag-cords: y v x");
                
            hold off 

        figure6 = figure;
            hold on
                xlabel("abs(x-cord)");
                ylabel("abs(y-cord)");
                zlabel("p-value");
                title("maximizing coords: abs");
                plot3(x, y, p);
        
            hold off 

        figure7 = figure;
            hold on
                xlabel("abs(x-cord)");
                ylabel("abs(y-cord)");
                zlabel("p-value");
                title("3D maximizing coords: abs");
                plot3(x, y, p);
            hold off

        figure8 = figure;
            hold on
                xlabel("Re(x-cord)");
                ylabel("Re(y-cord)");
                zlabel("p-value");
                title("3D maximizing coords: Re");
                plot3(xRe, yRe, p);
            hold off

        figure9 = figure;
            hold on
                xlabel("Im(x-cord)");
                ylabel("Im(y-cord)");
                zlabel("p-value");
                title("3D maximizing coords: Im");
                plot3(xIm, yIm, p);
            hold off

        figure10 = figure;
        hold on
            xlabel("Re(x-cord)");
            ylabel("Im(x-cord)");
            zlabel("p-value");
            title("3D maximizing coords: x Re vs Im");
            plot3(xRe, xIm, p);
        hold off

        figure11 = figure;
        hold on
            xlabel("Re(y-cord)");
            ylabel("Im(y-cord)");
            zlabel("p-value");
            title("3D maximizing coords: y Re vs Im");
            plot3(yRe, yIm, p);
        hold off
    end 
    
    fprintf('\nProgram Runtime:\n');
    toc(programTime)
    
%     saveas(figure2,'randomMatrix32_.1_10_20_hyg_1_10_.1.fig' ) 
end 