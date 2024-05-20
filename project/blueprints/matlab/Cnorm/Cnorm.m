classdef Cnorm
    properties
        
        modMatrix;
        cMatrix;
        hMatrix;
        
    end 

    methods(Static)

        function obj = Cnorm(N)
            obj.cMatrix = complex(rand(N, N), rand(N,N));
            obj.modMatrix = abs(obj.cMatrix);
            obj.hMatrix  = transpose(conj(obj.cMatrix)) * obj.cMatrix;

        end

    end

    methods(Static)

        function res = P1(obj)
            tic
            oneVector = ones(1,length(obj.modMatrix(1,:)));
            % modMatrix = sqrt((real(cMatrix)) .^ 2 + (imag(cMatrix)) .^ 2);
            % on avg saves a bit of time compared to above line
            normColumn =  oneVector * obj.modMatrix;
            [res, ~] = max(normColumn);
            toc
        end

        
        function res = P2power(obj,err_a)
            tic
            loc = obj.colMaxP1_P2(obj);
            vNow = obj.hMatrix(:, loc);
            
            error = 1;
            res = 0;
            % timeOutCounter = 0;
            
            while (error > err_a)

                vNext = obj.hMatrix * vNow;
        
                for i=1:length(obj.hMatrix(1,:))
                    if abs(vNow(i))
                        break;
                    end 
                end
                guess = vNext(i) / vNow(i);
                error = abs(res - guess);
                res = guess;
                vNow = vNext;
                % timeOutCounter = timeOutCounter + 1;
                % if (timeOutCounter > 10000)
                %     break;
                % end 
                

                
            end 
            
            res = sqrt(guess);
            toc


        end 

        function res = P1_P2(obj)
            tic
            modSquared = (obj.modMatrix) .^ 2;
            normColumn =  ones(1,length(obj.modMatrix(1,:))) * modSquared;
            [res, ~] = max(normColumn);
            if (length(res) - 1)
                res = res(1);
            end 

            res = sqrt(res);
            toc


        end 

        function res = colMaxP1_P2(obj)
            modSquared = (obj.modMatrix) .^ 2;
            normColumn =  ones(1,length(obj.modMatrix(1,:))) * modSquared;
            [~, res] = max(normColumn);
            if (length(res) - 1)
                res = res(1);
            end


        end 
 
        

        

    end


end 
