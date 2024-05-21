function pCompPlot(N, pMax)
    
    cMatrix = Cnorm(N);
    i = 1:pMax;
    norms = zeros(1, length(i));


    for j = 1:length(i)
        norms(j) = Cnorm.pPower(cMatrix, i(j), .0000001);
    end
    norms(j) = Cnorm.pPower(cMatrix, inf, .0000001);
    figure()
        
        xlabel("p");
        ylabel("Operator Norm Value");
        title("Operator Norm vs p");
        plot(i, norms);

end 