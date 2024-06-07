%Write a function that takes a set of p and throws in p and all the
%holder's conjugates
function p = generatePValues(pIn)
    p = pIn;
    for i = 1:length(pIn)
        q = pIn(i)/(pIn(i)-1);
        if ~ismember(q, p) && pIn(i) ~= 1
            p = sort([p, q]);
        end
    end
end


%Write a fucntion that returns the next index of a given value, increasing
%array
function index = findIndex(pIn, val)
    for t = 1:length(pIn)
        if pIn(t) < val
            t = t+1;
        elseif pIn(t) == val
            index = 0;
            return;
        else
            index = t;
            return;
        end
    end
end

        
%Figure out if there's a way to do this easier in matlab
%How to check if value is in an array in matlab?