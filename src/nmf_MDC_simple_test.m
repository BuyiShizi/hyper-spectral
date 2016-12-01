%% nmf_MDC_simple test script
    clear all
    % generate true W, H, V
    data_size = 100;
    W_true = abs( randn( data_size, 3 ) );
    sum_W = sum( W_true, 2 );
    for i = 1:data_size;
        W_true( i, : ) = W_true( i,: ) ./ sum_W(i);
    end
    H_true = abs( randn( 3, 2 ) );
    V = W_true * H_true;
    
    % factorize V using nmf_MDC_simple method
    [ W_test, H_test ] = nmf_MDC_simple( V, 3, 0.01, 10000 );
    
    M = [ 0 1 0; 0 0 1; 1 0 0];
    D_true = H_true - ( M * H_true);
    D_true = sum( sum( D_true.^2 ) );
    D_test = H_test - ( M * H_test);
    D_test = sum( sum( D_test.^2 ) );
    
    disp 'true distance',D_true
    disp 'test distance',D_test
    e_H = H_true - H_test;
    e_H_sum = sum( sum( e_H.^2 ) )
    V_test = W_test * H_test;
    e_V = V( :, 1:2 ) - V_test;
    e_V_sum = sum( sum( e_V.^2 ) )
    
    % visualize rest
    scatter( V(:,1), V(:,2), 'r' ); hold on ; scatter( V_test(:,1), V_test(:,2), 'k' );
    scatter( H_true(:,1), H_true(:,2) ,'filled','r')
    scatter( H_test(:,1), H_test(:,2) ,'filled','k')
    
%% nmf_MDC_constraint test
    clear all
    % generate true W, H, V
    data_size = 100;
    W_true = abs( randn( data_size, 3 ) );
    sum_W = sum( W_true, 2 );
    for i = 1:data_size;
        W_true( i, : ) = W_true( i,: ) ./ sum_W(i);
    end
    H_true = abs( randn( 3, 2 ) );
    V = W_true * H_true;
    
    % factorize V using nmf_MDC_simple method
    [ W_test, H_test ] = nmf_MDC_constraint( V, 3, 0.01, 0.8, 0.01, 30000 );
    
    M = [ 0 1 0; 0 0 1; 1 0 0];
    D_true = H_true - ( M * H_true);
    D_true = sum( sum( D_true.^2 ) );
    D_test = H_test - ( M * H_test);
    D_test = sum( sum( D_test.^2 ) );
    
    disp 'true distance',D_true
    disp 'test distance',D_test
    e_H = H_true - H_test;
    e_H_sum = sum( sum( e_H.^2 ) )
    V_test = W_test * H_test;
    e_V = V( :, 1:2 ) - V_test;
    e_V_sum = sum( sum( e_V.^2 ) )
    
    % visualize rest
    scatter( V(:,1), V(:,2), 'r' ); hold on ; scatter( V_test(:,1), V_test(:,2), 'k' );
    scatter( H_true(:,1), H_true(:,2) ,'filled','r')
    scatter( H_test(:,1), H_test(:,2) ,'filled','k')
    
%% nmf_MVC test
    clear all
    % generate true W, H, V
    data_size = 100;
    N = 3; % endmember number
    band_num = 2; % band number
    W_true = abs(randn(N, data_size));
    sum_W = sum(W_true);
    
    for i = 1:data_size
        W_true(:, i) = W_true(:, i) ./ sum_W(i);
    end
    H_true = abs(randn(band_num, N));
    V_true = H_true * W_true;
    
    % factorization V_true using nmf_MVC methos
    HasInit = 1;
    sigma = 0.5; % iteration tolerance
    tho = 0.8; % scaling factor that stepsize is reduced in every iteration
    alpha0 = 0.01; % initial step size for H
    beta0 = 0.01; % initial step size for W
    delta = 5; % positive number to control the effect of the sum-to-one
                 % constraint.
    lamda = 1; % control the tradeoff between the accurate reconstruction
                 % abd the volumn constraint
    tolerance = 0.1;
    iter_max = 10000;
    H_init_index = n_findr(V_true', 3);
    [W_test, H_test, H_record, Error] = nmf_MVC(V_true, N, H_init_index, HasInit, ...
                               sigma, tho, alpha0, beta0, ...
                               delta, lamda, ...
                               tolerance, iter_max);
    iter_max = length(Error);
    H_test = H_test(1:band_num, :);
    V_test = H_test * W_test;
    
    % visulize reconstruction
    figure_index = 1;
    figure(figure_index);
    figure_index = figure_index + 1;
    scatter(V_true(1, :), V_true(2, :), 'r'); hold on;
    scatter(V_test(1, :), V_test(2, :), 'k');
    scatter(H_true(1, :), H_true(2, :), 'filled', 'r');
    scatter(H_test(1, :), H_test(2, :), 'filled', 'k');
    plot(reshape(H_record(1, 1, :), 1, iter_max), ...
         reshape(H_record(2, 1, :), 1, iter_max));
    plot(reshape(H_record(1, 2, :), 1, iter_max), ...
         reshape(H_record(2, 2, :), 1, iter_max));
    plot(reshape(H_record(1, 3, :), 1, iter_max), ...
         reshape(H_record(2, 3, :), 1, iter_max));
    
    figure(figure_index)
    figure_index = figure_index + 1;
    plot(Error);
