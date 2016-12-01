%% function to solve least linear square
%-----------------------------------------------------------
    %this function can solve constrained linear equations
    %with constrains of the sum of coefficient is 1;
    %the solution is that we add the equation of sum to 1 
    %to the linear equation.
%-----------------------------------------------------------
function coefficient = linear_least_square (X,y)
    % X:            input metrix
    % y:            input vector
    % coefficient:  output solution
    
    % we need to add the constrains of sum of coefficient is 1
    % the constrained linear equation became unconstrained linear equation:
    %
    % | 1   1   1   1   1 |   | 1 |
    % | x   x   x   x   x |   | y |
    % | x   x   x   x   x | = | y |
    % | x   x   x   x   x |   | y |
    % | x   x   x   x   x |   | y |
    % | x   x   x   x   x |   | y |
    [row_X, col_X] = size (X);
    X = [ X];
    y = [y];
    
    % now we calculate the solutions
    X_square_metrix = X' * X;
    X_square_inv = inv (X_square_metrix);
    
    coefficient = X_square_inv * X' * y;
end

%% how to solver endmember extraction by using this function
% when solve endmember extraction:
% X represents endmember data, 
% the row of X represents all endmember data in one band;
% the column of X represents all bands of data of one endmember;
% y is a vector represents all band data of one pixel;
% so when we do endmember extraction with this function, we need traverse all the pixels in hyperspectral image.
