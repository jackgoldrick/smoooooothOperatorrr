classdef Cnorm
    properties
        
        modMatrix;
        
    end 

    methods(Static)

        function obj = Cnorm(N)
            cMatrix = complex(rand(N, N), rand(N,N));
            obj.modMatrix = abs(cMatrix);

        end

    end

    methods(Static)

        function res = P1(obj)
            tic
            oneVector = ones(1,N);
            % modMatrix = sqrt((real(cMatrix)) .^ 2 + (imag(cMatrix)) .^ 2);
            % on avg saves a bit of time compared to above line
            normColumn =  oneVector * obj.modMatrix;
            [res, ~] = max(normColumn);
            toc
        end

        
        function res = P2(obj)
        

        end 

        function res = P1_P2(obj)
            tic
            modSquared = (obj.modMatrix) .^ 2;
            normColumn =  oneVector * modSquared;
            [res, ~] = max(normColumn);
            if (length(res) - 1)
                res = res(1);
            end 

            res = sqrt(res);
            toc


        end 

        

    end


end 
