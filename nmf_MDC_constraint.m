function [ W, H ] = nmf_MDC_constraint( V, N, lamda0, lamda1, Tolerance, IterMax )
%% Minimum distant constraint for 
%  nonnegative matrix factorization
% W,H:       output for factorized matrixs
% V:         original lined hyper data with some bands
% N:        endmember number
% alpha0:    distance constraint factor
% alpha1:    sum to 1 constraint factor
% Tolerance: max error when to stop iteration
% InterMax:  max interation number for stop iteration

%% start process

    % V( row_V, col_V ) = W_cur( row_V, N) * H_cur( N, col_V )
    P = [ 0 1 1; 1 0 1; 1 1 0 ];
    I = [ 1 0 0; 0 1 0; 0 0 1 ];
    [ row_V, col_V ] = size( V );
    V = [ V, lamda1*ones(row_V, 1) ];
    W_cur = zeros( row_V, N );
    H_cur = zeros( N, col_V );

    % initial W, H
    W = abs( randn( row_V, N ) );
    sum_W = sum( W, 2 );
    for i = 1:row_V
        W( i, : ) = W( i,: ) ./ sum_W(i);
    end
    
    H = abs( randn( N, col_V ) );
    H = [ H lamda1*ones(N,1) ];
    
    % calculate linear least square error
    e = (V - (W * H)) ./ V;
    e_sum = sum( sum( e.^2 ) );
    iter_cur = 1;
    
    % start iteration
    while( e_sum > Tolerance )
        if( iter_cur == IterMax )
            break;
        end
        
        % update H
        normi = (W' * V) + (lamda0 * P * H);
        denormi = ( (W' * W) + (lamda0 * I) ) * H;
        update_factor = normi ./ denormi;
        H = H .* update_factor;
        H( :,col_V+1 ) = lamda1 * ones( N, 1 );
        
        % update W
        normi = V * H';
        denormi = W * H * H';
        update_factor = normi ./ denormi;
        W = W .* update_factor;
        
        % calculate error
        e = (V - ( W * H )) ./ V;
        e_sum = sum( sum( e.^2 ) );
        
        iter_cur = iter_cur + 1;
    end
    
    H = H( :,1:2 );
end

