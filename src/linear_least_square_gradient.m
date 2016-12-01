%% This is a function for linear least square problem solved by
%   gradient descent method
function [X, E] = linear_least_square_gradient(Y, C, X_INIT, RATE, TOLERANT, TOLERANT1)
% X:                    output for abundance solution
% E:                    output final solution error
% Y:                    input for object value
% C:                    input for abundance data
% X_INIT:               input for initial value of abundance solution
% RATE:                 update rate in gradient descent
% TOLERANT:             tolerance to control when to stop iteration

	variable_num = length(X_INIT); % get the number of endmember
    y_diff_x = ones (variable_num,1); % store gradient
    max_num = 1e4;
    
    %% plot error figure
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
    sum_Y = max(Y);
    Y = Y/sum_Y;
    C = C/sum_Y;
    %% 
    x = X_INIT;
    rate = RATE;
    %y_diff_x = 2 * x * C * C' - 2 * Y * C';
    y_diff_x = 2 * (C' * C) * x - 2 * C' * Y;
    max_diff = sum(abs(y_diff_x));
    diff_refer_orig = max_diff;
    diff_refer = diff_refer_orig;
    
    f_error = Y - C * x;
    e_previous = 0;
    e = sum (f_error .^2 );
    i = 1;
    
    diff_rate = 1;
    while (e > TOLERANT)
%         y_diff_x = 2 * x * C * C' - 2 * Y * C';
        y_diff_x = 2 * (C' * C) * x - 2 * C' * Y;
        x_pre = x;
        max_diff = max(abs(y_diff_x));
        x = x - rate * y_diff_x % update abundance
        y_diff_x = 2 * (C' * C) * x - 2 * C' * Y;
%         quiver(x(1), x(2), y_diff_x(1)/max_diff, y_diff_x(2)/max_diff);
%         y_diff_x = 2 * x * C * C' - 2 * Y * C';
        
        f_error = Y - C * x; % error for current abundance
        
        e_previous = e;
        e = sum ( f_error.^2 )
        max_diff_pre = max_diff;
        max_diff = sum(abs(y_diff_x))
        diff_rate_pre = diff_rate;
        diff_rate = abs(e-e_previous)/e_previous;
        if(abs(diff_rate-diff_rate_pre)<1e-3)
            rate = rate * 2;
        end
        if (abs(diff_rate-diff_rate_pre) > 0.2)
            rate = rate /2;
        end
        
        % change rate
%         i = 1;
        if (max_diff < (diff_refer/10))
            i = i+1;
%             rate = rate/(1);
            rate = rate;
            diff_refer = diff_refer/10;
%             i = i + 1;
        end
    end
    
    X = x;
    E = e;  
end % end function
