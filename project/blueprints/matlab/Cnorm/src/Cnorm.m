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
        %% P1 -> P1
        function res = P1(obj)
            tic
            oneVector = ones(1,length(obj.modMatrix(1,:)));
            % modMatrix = sqrt((real(cMatrix)) .^ 2 + (imag(cMatrix)) .^ 2);
            % on avg saves a bit of time compared to above line
            normColumn =  oneVector * obj.modMatrix;
            [res, ~] = max(normColumn);
            toc
        end
        %% Pinf -> Pinf
        function res = Pinf(obj)
            tic
            oneVector = ones(length(obj.modMatrix(1,:)), 1);
            % modMatrix = sqrt((real(cMatrix)) .^ 2 + (imag(cMatrix)) .^ 2);
            % on avg saves a bit of time compared to above line
            normColumn =  obj.modMatrix * oneVector;
            [res, ~] = max(normColumn);
            toc


        end 
        
        function res = pPower(obj, p, err_a)
            %% put this back to 2 when done testing
            if p == 0

                res = obj.P2power(obj, err_a);
                return

            elseif p == 1

                res = obj.P1(obj);
                return

            elseif p == inf

                res = obj.Pinf(obj);
                return

            end

            cMatrix = obj.cMatrix;
            cMatrixPrime = obj.cMatrix';

            tic
            loc = obj.colMaxP(obj, p);
            vNow = obj.hMatrix(:, loc);
            q = 1 / (1 - 1/p);
            error = 1;
            res = 0;
            oldGuess = Cnorm.vectorPNorm(vNow, p);
            vNow = vNow ./ oldGuess;
            while (error > err_a)

                vNext = cMatrix * vNow;

                guess = Cnorm.vectorPNorm(vNext, p);

                vNextDualNormed =  vNext ./ guess;

                z = cMatrixPrime * vNextDualNormed;

                error = abs((guess - oldGuess) / guess);

                % if Cnorm.vectorPNorm(z,q) <= z' * vNow| | error < err_a
                if error < err_a
                    break;

                end

                vNow =  z ./ Cnorm.vectorPNorm(z, q);
                oldGuess = guess;


                
            end 
            
            res = guess;
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
        
        function res = colMaxP(obj, p)
            modSquared = (obj.modMatrix) .^ p;
            normColumn =  ones(1,length(obj.modMatrix(1,:))) * modSquared;
            [~, res] = max(normColumn);
            if (length(res) - 1)
                res = res(1);
            end


        end

        function res = vectorPNorm(v, p)
            if (length(v(1,:)) == 1)
                res = ones(1,length(v)) * (abs(v) .^ p);
            else  
                res = (abs(v) .^ p) * ones(length(v), 1);
            end
            res = (res) ^ (1/p);
                
        end 
               

 
        

        

    end


end 
