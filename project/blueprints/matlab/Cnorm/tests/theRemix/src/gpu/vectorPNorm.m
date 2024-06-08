
%{
    This Function takes the p norm along dimension 1 and preserves the
    2nd and 3rd dimensions
  
%}

function res = vectorPNorm(v, p)
        if isinf(p) || isnan(p)
        

            res = max(abs(v),[], 1);
            return
        end

        res = sum(abs(v) .^ p, 1);
        res = (res) .^ (1/p);
        % if row Vector
        % if (length(v(:,1)) == 1)
        %     res = (abs(v) .^ p) * ones(length(v(1,:)), 1);
        % else  
        %     res = ones(1,length(v(:,1))) * (abs(v) .^ p);
        % end
        % res = (res) .^ (1/p);


            
    end