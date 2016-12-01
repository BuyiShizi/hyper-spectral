function [W, H, H1, H2, H3, E, IT] = hyper_nmf (V, N, alpha, tolerance, maxIter)
%% this function use nmf-similiar method to factorization hyperspectrum
% W,H:          output----factored output for hyperspectral data
% V:            input ----original hyperspectroal
% N:            input ----endmember numbers
% alpha:        input ----scale factor 
% tolerance:    input ----input for error tolerance
% MaxIter:      input ----input for max iteration times

%%

[row_V, col_V] = size ( V );
V_EXT = [alpha*V, ones(row_V,1)];
W = rand ( row_V, N );
sum_W = sum (W, 2);
for i = 1:row_V;
    W(i,:) = W(i,:) ./ sum_W(i);
end
% Hi = randi ( row_V );
H = [zeros(N, col_V),ones(N,1)]; % add sum is 1 to the equation
for i = 1:N
    Hi = randi ( row_V );
    H(i,:) = V_EXT(Hi,:);
end

H1 = zeros(2,row_V);
H2 = zeros(2,row_V);
H3 = zeros(2,row_V);

% record error for every iteration
E = zeros(1,maxIter);
% get the error

e = V_EXT - W * H;
e2 = sum ( sum (e .^2) );
Iter_n = 1; % iteration times
j = 0; % recorde iter

    while ( e2 > tolerance )
        if (Iter_n == maxIter )
            break;
        end

        % get nex H
        WH = W * H;
        V_WH = V_EXT ./ WH;
        [r_, c_] = size(H);
        for i_ = 1:r_
            for j_ = 1:c_
                WV_WH = sum( W(:,i_).*V_WH(:,j_) );
                Wi_sum = sum( W(:,i_) );
                H(i_,j_) = H(i_,j_) * WV_WH / Wi_sum;
            end
        end

        % get nex W
        [r_, c_] = size(W);
        for i_ = 1:r_
            for j_ = 1:c_
                HV_WH = sum( H(j_,:).*V_WH(i_,:) );
                Hi_sum = sum( H(j_,:) );
                W(i_,j_) = W(i_,j_) * HV_WH / Hi_sum;
            end
        end
        

        % get error
        e = V_EXT - W * H;
        e2 = sum ( sum (e .^2) );

        % recorde iteration
        
        j = j + 1;
        H1(1,j) = H(1,1);
        H1(2,j) = H(1,2);
        H2(1,j) = H(2,1);
        H2(2,j) = H(2,2);
        H3(1,j) = H(3,1);
        H3(2,j) = H(3,2);
        
        E(Iter_n) = e2;
            
        Iter_n = Iter_n + 1;
    end
    
    H1 = H1(:,1:j);
    H2 = H2(:,1:j);
    H3 = H3(:,1:j);
    E = E(1:Iter_n);
    IT = Iter_n;
end