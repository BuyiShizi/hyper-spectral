%% This is a function for linear least square problem solved by
%   gradient descent method
function [X, E] = gradient_descent(Y, C, X_INIT, RATE, TOLERANT)
% X:                    output for abundance solution
% E:                    output final solution error
% Y:                    input for object value
% C:                    input for abundance data
% X_INIT:               input for initial value of abundance solution
% RATE:                 update rate in gradient descent
% TOLERANT:             tolerance to control when to stop iteration

	variable_num = length(X_INIT); % get the number of endmember
    y_diff_x = ones (variable_num,1); % store gradient
    
    %% plot error figure for gradient observation
%     figure(1)
%     hold on
%     [x1,x2] = meshgrid (0:0.01:1,0:0.01:1);
%     [x1_size, x2_size] = size(x1);
%     e = zeros (x1_size, x2_size);
%     for i = 1:x1_size
%         for j = 1:x2_size
%             f_error = Y - C * [x1(i,j); x2(i,j)];
%             e(i,j) = sum (f_error.^2);
%         end
%     end
%     contour(x1, x2, e,50); 
    %% renormalize input data
    % normalize the input data is to reduce the range of gradient
%     sum_Y = max(Y);
%     Y = Y/sum_Y;
%     C = C/sum_Y;
    %% 
    x = X_INIT;
    rate = RATE;
    y_diff_x = 2 * (C' * C) * x - 2 * C' * Y; 
    
    %% calculate error for current abundance
    f_error = Y - C * x;
    e_previous = 0;
    e = max(sum (f_error .^2 ));  
    diff_rate = 1; % used for chage step size
    while (e > TOLERANT)
        y_diff_x = 2 * (C' * C) * x - 2 * C' * Y;
        x = x - rate * y_diff_x % update abundance
        y_diff_x = 2 * (C' * C) * x - 2 * C' * Y;
%         quiver(x(1), x(2), y_diff_x(1)/max_diff, y_diff_x(2)/max_diff);
        
        f_error = Y - C * x; % error for current abundance
        
        e_previous = e;
        e = max(sum ( f_error.^2 ))
        diff_rate_pre = diff_rate;
        diff_rate = abs(e-e_previous)/e_previous;
        
        % if rate is too small, enlarge it
        if(abs(diff_rate-diff_rate_pre)<1e-3)
            rate = rate * 3;
        end
        % if rate is too big to deverge, small it
        if (abs(diff_rate-diff_rate_pre) > 0.2)
            rate = rate /3;
        end
        
    end
    
    X = x;
    E = e;  
end % end function

