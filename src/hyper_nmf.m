function [W, H, H1, H2, H3, E, IT] = hyper_nmf (V, N, H_I, alpha, tolerance, maxIter)
%% this function use nmf-similiar method to factorization hyperspectrum
% W,H:          output----factored output for hyperspectral data
% V:            input ----original hyperspectroal
% N:            input ----endmember numbers
% H_I:			H initial value
% alpha:        input ----scale factor 
% tolerance:    input ----input for error tolerance
% MaxIter:      input ----input for max iteration times

%%
[row_V, col_V] = size ( V );
V_EXT = [alpha*V, ones(row_V,1)];

V_cut = V_EXT(1:1:row_V,:);
[row_V, col_V] = size ( V_cut );
W = abs( randn ( row_V, N ));
sum_W = sum (W, 2);
for i = 1:row_V;
    W(i,:) = W(i,:) ./ sum_W(i);
end
% Hi = randi ( row_V );
H = [H_I,ones(N,1)]; % add sum is 1 to the equation
% for i = 1:N
%     Hi = randi ( row_V );
%     H(i,:) = V_cut(Hi,:);
% end

H1 = zeros(3,row_V);
H2 = zeros(3,row_V);
H3 = zeros(3,row_V);

j = 1; % recorde iter
H1(1,j) = H(1,1)*alpha;
H1(2,j) = H(1,2)*alpha;
H1(3,j) = H(1,3)*alpha;
H2(1,j) = H(2,1)*alpha;
H2(2,j) = H(2,2)*alpha;
H2(3,j) = H(2,3)*alpha;
H3(1,j) = H(3,1)*alpha;
H3(2,j) = H(3,2)*alpha;
H3(3,j) = H(3,3)*alpha;

% record error for every iteration
E = zeros(1,maxIter);
% get the error

e = V_cut - W * H;
e2 = sum ( sum (e .^2) ) ;
Iter_n = 1; % iteration times

    while ( e2 > tolerance )
        if (Iter_n == maxIter )
            break;
        end

        % get nex H
        % keep H for one iteration
        H_k = H;
        WV = W' * V_cut;
        WWH = W' * W * H;
        H = H .* ( WV ./ WWH ); % update H

        % get nex W
        VH = V_cut * H';
        WHH = W * ( H * H' );
        W = W .* ( VH ./ WHH ); % update W
        max_W = max(W);
        i1_ = find( W(:,1) == max_W(1) );
        i2_ = find( W(:,2) == max_W(2) );
        i3_ = find( W(:,3) == max_W(3) );
        W(i1_,1) = 1;
        W(i2_,2) = 1;
        W(i3_,3) = 1;
        W_sum = sum(W, 2);
        [r_, c_] = size(W);
        for i_ = 1:c_
            W(:,i_) = W(:,i_) ./ W_sum;
        end

        % get error
        e = V_cut - W * H;
        e2 = sum ( sum (e .^2) );

        % recorde iteration
        
        j = j + 1;
        H1(1,j) = H(1,1);
        H1(2,j) = H(1,2);
        H1(3,j) = H(1,3);
        H2(1,j) = H(2,1);
        H2(2,j) = H(2,2);
        H2(3,j) = H(2,3);
        H3(1,j) = H(3,1);
        H3(2,j) = H(3,2);
        H3(3,j) = H(3,3);
        
        E(Iter_n) = e2;
            
        Iter_n = Iter_n + 1;
    end
    
    H1 = H1(:,1:j);
    H2 = H2(:,1:j);
    H3 = H3(:,1:j);
    E = E(1:Iter_n);
    IT = Iter_n;
end