function genCompPlot(N, dp, pMax, type)

    switch type
        case 'h'
            cMatrix = hadamard(N);
        case 'r'
            cMatrix = rand(N, N);
        case 'c'
            cMatrix = complex(rand(N, N), rand(N,N));
        case 'a'
            cMatrix = [ 0 1 -1; -1 0 1; 1 -1 0];
            cMatrix = cMatrix .* (1 / sqrt(3));
        case 'v'
            cMatrix = [ 1 0 1 0; 0 1 0 1];
    end

    tic


    p = 1:dp:pMax;

    norms = zeros(1, length(p));
    correctNorms =zeros(1, length(p));


    for j = 1:length(p)
        q = 1 / (1 - 1 / p(j));
        
        norms(j) = genComparison(cMatrix, p(j), .0000001);

        if type == 'h'
            correctNorms(j) = max(N ^ (1/p(j)), N ^ (1 / q));
        end 

        if type == 'v'
            correctNorms(j) = 2 ^ (1 / q);
        end

        if p(j) <= 2 && type == 'a'
            correctNorms(j) = ((1 + 2 ^(p(j) - 1)) ^ (1/(p(j)))) / sqrt(3);

        elseif p(j) > 2 && type == 'a'
            correctNorms(j) = ((1 + 2 ^(q - 1)) ^ (1/(q))) / sqrt(3);
        end 


    end

    % norms(j) = hadComparison(cMatrix, inf, .0000001);
    figure
        hold on
        xlabel("p");
        ylabel("Operator Norm Value");
        title("Operator Norm vs p");
        plot(p, norms, '--b', 'LineWidth', 1);
        plot(p, correctNorms, '-r');
        legend("Our Method", "Exact Value");
        hold off
    toc

end 