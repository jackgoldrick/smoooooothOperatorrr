function averageProportion = plotBMaxEffectiveness(N, dp, pMax, sMax, invertableMatrix, DSBMaxesDesired, trials)
    proportions = plotMaxSimDiag(N, dp, pMax, sMax, invertableMatrix, DSBMaxesDesired);
    for j = 2:trials
        newProportions = plotMaxSimDiag(N, dp, pMax, sMax, invertableMatrix, DSBMaxesDesired);
        proportions = cat(2, proportions, newProportions);
    end
    averageProportion = mean(proportions);

    edges = [0, .8, 1.2, 2, 4, 8];

    figure;
    histogram(proportions, edges);
    
    xlabel('proportion');
    ylabel('Frequency');
end

    