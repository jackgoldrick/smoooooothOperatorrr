    function res = vectorPNorm(v, p)
        if isinf(p)
            res = max(abs(v));
            return
        end
        
        if (length(v(1,:)) == 1)
            res = ones(1,length(v)) * (abs(v) .^ p);
        else  
            res = (abs(v) .^ p) * ones(length(v), 1);
        end
        res = (res) .^ (1/p);
            
    end