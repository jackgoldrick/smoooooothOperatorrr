function pDimCompPlot(N, ~)
    cMatrix = Cnorm(N);
    i = 1:pMax;
    norms = zeros(1, length(i));


    for j = 1:length(i)
        norms(j) = Cnorm.pPower(cMatrix, i(j), .0000001);
    end

    figure()
        hold on
            xlabel("Dimension");
            ylabel("Operator Norm Value");
            title("Operator Norm vs Size");
            plot(i, norms);
        hold off

end 














