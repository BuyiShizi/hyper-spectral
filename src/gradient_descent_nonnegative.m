function [X, E] = gradient_descent_nonnegative(Y, C, MaxIter, Tolerance)
%% thid functio solve nonnegative linear least square method
%  

%% do some praparation 
[row_Y, col_Y] = size(Y)
[row_C, col_C] = size(C)
X = rand(col_C, col_Y); % initialize X with random variables
%% calculate first iteration
iter_times = 1;
error_matrix = Y - C * X;
error_matrix = error_matrix .^ 2;
E = sum(sum(error_matrix));

%% start iteration
while(E > Tolerance)
    if (iter_times == MaxIter)
        break;
    end
    
    v1 = C' * Y;
    v2 = C' * C * X;
    scale_factor = v1 ./ v2;
    
    % update abundance using scale_factor
    X = X .* scale_factor;
    
    error_matrix = Y - C * X;
    error_matrix = error_matrix .^ 2;
    E = sum(sum(error_matrix));
    iter_times = iter_times + 1;
end

end