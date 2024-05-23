function genCompPlot(N, dp, pMax, type, g)
    %% Add you matrix you want to test here
    switch type
        case 'h'
            cMatrix = hadamard(N);
            norm(cMatrix, 1)
        case 'r'
            if (length(N) - 1)
                cMatrix = rand(N(1), N(2));

            else
                   cMatrix = rand(N, N);
            end
        case 'c'
            
             if (length(N) - 1)
                  cMatrix = complex(rand(N(1), N(2)), rand(N(1),N(2)));
            else
                 cMatrix = complex(rand(N, N), rand(N,N));
            end
        case 'm'
            cMatrix = [ 0 1 -1; -1 0 1; 1 -1 0];
            cMatrix = cMatrix .* (1 / sqrt(3));
        case 'a'
            cMatrix = [ 1 0 1 0; 0 1 0 1];
        case 'p'
            cMatrix = transpose([ 1 0 1 0; 0 1 0 1]);
        case 'd'
            cMatrix = [ -2 1 1; 1 -2 1; 1 1 -2];
            cMatrix = cMatrix .* (1 / -3);
    end

    tic


    p = 1:dp:pMax;

    norms = zeros(1, length(p));
    if type == 'a' || type == 'p'
        correctNormsMax = zeros(1, length(p));
        correctNormsQ = zeros(1, length(p));
        correctNormsP = zeros(1, length(p));

    else 
        correctNorms = zeros(1, length(p));
    end
    % Q = zeros(1, length(p));

    %% Add  your solution here  
    for j = 1:length(p)
        q = 1 / (1 - 1 / p(j));
        % Q(j) = q;
        
        norms(j) = genComparison(cMatrix, p(j), .0000001);

        if type == 'h'
            correctNorms(j) = max(N ^ (1/p(j)), N ^ (1 / q));
     

        elseif type == 'a' || type == 'p'
            correctNormsMax(j) = max(2 ^ (1 / q), 2 ^ (1 / p(j)));
            correctNormsQ(j) = 2 ^ (1 / q);
            correctNormsP(j) = 2 ^ (1 / p(j));

        elseif p(j) <= 2 && type == 'm'
            correctNorms(j) = ((1 + 2 ^(p(j) - 1)) ^ (1/(p(j)))) / sqrt(3);

        elseif p(j) > 2 && type == 'm'
            correctNorms(j) = ((1 + 2 ^(q - 1)) ^ (1/(q))) / sqrt(3);
        end 

         


    end

    if type == 'c'|| type == 'r'
        
        for i=2:length(norms)
            if ~(isnan(norms(i)))
                break;
            end
        end
        % I wonder if this ratio can be seen as an approximation error
        % need to compute sequence of ratios.
        % this should converge 
        if i == length(norms)
            ratio = 0 %#ok<NASGU,NOPRT>
        else
            i %#ok<NOPRT>
            ratio = norms(i) / norms(1) %#ok<NASGU,NOPRT>
            pValue = 1 + i *dp %#ok<NASGU,NOPRT>

        end

            
        
        
    end
    % norms(j) = hadComparison(cMatrix, inf, .0000001);

    if g == 'y'
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
                    
                elseif type == 'd'
                    plot(p, norms, '--b', 'LineWidth', 1);
                    plot(p, 1, '-r');
                    legend("Our Method", "Lower Bound");
                else 
                    plot(p, norms, '--b', 'LineWidth', 1);
                    plot(p, correctNorms, '-r');
                    legend("Our Method", "Exact Value");
                end
            hold off
    end
    toc

end 