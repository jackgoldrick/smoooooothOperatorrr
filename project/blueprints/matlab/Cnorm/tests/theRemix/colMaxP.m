        function res = colMaxP(modMatrix, p)
            modSquared = (modMatrix) .^ p;
            normColumn =  ones(1,length(modMatrix(:,1))) * modSquared;
            [~, res] = max(normColumn);
            if (length(res) - 1)
                res = res(1);
            end


        end