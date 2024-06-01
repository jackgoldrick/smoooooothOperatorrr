    function res = vectorPNorm(v, p)
        if isinf(p)
            res = max(abs(v));
            return
        end
        % if row Vector
        if (length(v(:,1)) == 1)
            res = (abs(v) .^ p) * ones(length(v(:,1)), 1);
        else  
            res = ones(1,length(v(:,1))) * (abs(v) .^ p);
        end
        res = (res) .^ (1/p);
            
    end