%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%
    %%                                                                                                                         %%
        %% Author: @jackgoldrick @wiwu2390                                                                                 %%
        %% Repository: https://github.com/jackgoldrick/smoooooothOperatorrr                                                %%
        %%                                                                                                                 %%                                                                                            
        %% 'y' - row or column vector that we want to find the p normalized dual                                           %%
        %% 'p' - operator norm                                                                                             %%
    %%                                                                                                                         %%
%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%

function yDual = dual(y, p)
    if p == 1
        yDual = conj(sign(y));
        return
    elseif isinf(p) || isnan(p)
        yDual = 0 .* y;
        [val,ind] = max(y, [], 1);
        loc = size(y);
        yDual(1:loc(1) == ind) = conj(sign(val));
        return
    end 
    q = 1 / (1 - 1/p);

    yDual = (abs(y).^ (p-1)) .* (sign(y));
    
    yDual =  yDual ./ vectorPNorm(yDual, q);
    
end 