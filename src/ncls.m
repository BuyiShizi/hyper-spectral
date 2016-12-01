%% linear least square with great than 0 constraints
function abundance_solution = ncls (A, y, N, tolerance)
    % A:                    input for coefficient metrix.
    % y:                    input for linear equation output.
    % N:                    input for number of variables for abundance_solution.
    % tolerance:            input for tolerance to stop iteration
    % abundance_solution:   output for solution for linear equation.
    
    %% initilize step
    
    % passive and active set
    P = zeros (1,N);
    R = 1:N;
    
    % initilize abundance solution
    d = zeros (N, 1);
    s = zeros (N, 1);
    
    % intermidiate variable
    w = A' * (y - A*d);
    
    %% main loop
    while ((~isempty(R)))
        w_max = 0;
        m = 0;
        [w_max, m] = max (w(R));
        Ri =  R(m);
        
        % include the index m in P and remove it from R
        P_prev = P;
        P(R(m)) = R(m);
        R(m) = 0;
        R = R(R~=0);
        
        % calculate intermidiate variable s
        A_P = A(:,P~=0);
        s_P = (A_P' * A_P) \ (A_P' * y);
        %s = zeros (N,1);
        s(P~=0) = s_P;
        %% inner loop
        while (min(s_P) < -tolerance) 
            tmp = d(P~=0) ./ (d(P~=0) - s_P);
            a = - min(tmp);
            d = d + a * (s - d);
            P = zeros (N,1);
            P(d>0) = find (d > 0);
            R = find (d <= 0);
            A_P = A (:,P~=0);
            s_P = (A_P' * A_P) \ (A_P' * y);
            s = zeros (N,1);
            s(P~=0) = s_P;
        end
        d = s;
        w = A' * (y - A * d);
    end
    
    %% return abundance solution
    abundance_solution = d; 
end
