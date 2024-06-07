function index = findIndex(pIn, val)
    for t = 1:length(pIn)
        if pIn(t) < val
            t = t+1;
        else
            index = t;
            return;
        end
    end
    index = 0;
end