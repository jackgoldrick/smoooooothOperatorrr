function plotMaxSimDiag(N, dp, pMax, sMax, invertableMatrix, DSBMaxesDesired, bound)

    if isempty(N) %#ok<BDLOG>
        cMatrix = [2 1; 1 2];
        N = 2;
    elseif (length(N) - 1) %#ok<BDLOG>
        cMatrix = complex(randn(N(1), N(2)), randn(N(1), N(2)));
    else 
         cMatrix = complex(randn(N, N), randn(N, N));
    end
           

    programTime = tic;
    p = 1:dp:pMax;
    sizeP = length(p);
 
    norms = zeros(sizeP, 1);

    %% Computes the P->P and/or P->R depending on pq flag
     % 'g' - computes P->R with O(sizeP * (sizeR + 1)) total operator norms calculated with graph
     % 'q' - computes P->R for ONE value of R ONLY
     % 'g' - computes P->P and P->R with O(sizeP * sizeR) total operator norms calculated without graph
     % 'p' - computes P->P ONLY 

    vMax2 = zeros(length(cMatrix(1,:)), 1);

    correctNorms = zeros(sizeP);
    
    if invertableMatrix == 0
        if (length(N)-1)
            N = N(2)
        end
        invertableMatrix = hadamard(N) ./ sqrt(N);
    end   

    %Needs to calculate once for first DSBMaxes set
    [norms(1), correctNorms(1), DSBMaxes] = maxSimDiag(cMatrix, invertableMatrix , p(j), 1e-1, sMax, 1e-7);
    
    for j = 2:sizeP
        [norms(j), correctNorms(j), DSBMaxes] = maxSimDiag(cMatrix, invertableMatrix , p(j), 1e-1, sMax + min(DSBMaxesDesired, 0), 1e-7, abs(DSBMaxesDesired), DSBMaxes);
    end


            
       

            
                for j = 1:sizeP

                    [norms(j), correctNorms(j), DSBMaxes] = maxSimDiag(cMatrix, hadamard(N) ./ sqrt(N), p(j), 1e-1, sMax, 1e-7, BMaxesDesired, DSBMaxes);
                    if bound == 'y'
                        normsBound1(j) = correctNorms(1)^(1/p(j)) * infNorm^(1-1/p(j));
                        if p(j) < 2
                            %We know these need to get fixed... our goal rn
                            %is to ignore though and come back to it
                            normsBound2(j) = correctNorms(1)^(2/p(j)-1) * twoNorm^(2-2/p(j));
                        else
                            normsBound2(j) = infNorm^(1-2/p(j)) * twoNorm^(2/p(j));
                        end        
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