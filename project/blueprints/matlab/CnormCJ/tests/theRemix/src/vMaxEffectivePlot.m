
%Goal: determine the relationship between dp size and accuracy for p, n.
%Plot the maximum dp size that yields norm within error tolerance % of the
%time, given n and p

%Plot largest dp possible against a and n

function vMaxEffectivePlot(N, pMax, dp, sMax, randInclusions, tolerance, accuracy)
%err=.0000001;N=[5,7];pMax=50;dp=.1;sMax=50;randInclusions=1;tolerance=.05;tests=100;
%accuracy=.1;
    
    err = .0000001;

    if length(N) ~= 2
        fprintf('Please Enter 2 values for N\n');
        return
    end
    
    %Creates n's of complex numbers to try
    n = N(1):N(2);

    %Creates p's to try for each n
    p = 1:dp:pMax;

    %Creates results recording - find max dp size yielding correct norm at
    %each tolerance using guessing method specified by previous vMax,
    %randInclusions #
    res = zeros(length(n), length(p));

    for j = 1:length(n)

        a = complex(rand(n(j),n(j)),rand(n(j),n(j)));
        
        minNorms = zeros(length(p),1);
        
        %Test value to make sure minNorms is being called correctly.
        vMax = zeros(n(j),length(p);

        %find the minNorms values here using the genComparison function
        %with the sMax  argument

        
        for k = 1:length(p)
            
            [minNorms(k), vMax(:,k)] = genComparison(a, p(k), err, sMax, 0);

            plot(p, minNorms, '-b', 'LineWidth', 1);


            %Find which dp works in this case that yields the best result
            %most of the time. Start with a given dp (brute force it in for
            %now) and write the code to check if it's oka, or no.
        end

        for k = 1:length(p)
            if k >= 2
                %p0 sets the previous p that you're capturing the
                %maximizing norm for
                count = 0;
                p0 = p(k-1);
                    
                [~, vTest] = genComparison(a, p0, err, sMax, 0);
                [nTest, ~] = genComparison(a, p(k), err, randInclusions, vMax);

                if(minNorms(k) - nTest > err)
                    count = count+1;
                end
            end
        end 

        