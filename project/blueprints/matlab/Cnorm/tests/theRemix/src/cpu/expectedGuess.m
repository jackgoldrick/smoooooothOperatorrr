function res = expectedGuess(x, y, z)
    res = 0;    
    for n = 1:min(z,y)
            res = res + nchoosek(y, n) * nchoosek(x-y, z-n);
    end
    res = res / nchoosek(x, z);
end


