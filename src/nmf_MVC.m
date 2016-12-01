function [W, H, H_record, E] = nmf_MVC(V, N, H_InitIndx, HasInit, ...
                          sigma, tho, alpha0, beta0, ... 
                          delta, lamda, ...
                          Tolerance, IterMax)
%% Minimum volumn constraint for 
%  nonnegative matrix factorization
% W,H:       output for factorized matrixs
% V:         original lined hyper data with some bands
% H_record:  record H in every iteration
% E:         record error in every iteration
% N:         endmember number
% H_InitIndx:init index for H according from nfindr
% HasInit:   flag indicate use H_InitIndx or not
% sigma:     iteration tolerance in MVC,[0, 1/2]
% tho:       scaling factor the the stepsize is reduced in every iteration
%            [0, 1]
% alpha0:    initial stepsize for update H
% beta0:     initial stepsize for update W
% delta:     factor to control the effect of sum-to-one constraint
% lamda:     factor to control the tradeoff between
%            the accurate reconstruction and the volumn constraint
% Tolerance: max error when to stop iteration
% IterMax:  max interation number for stop iteration

%% start process
  % data pre-process
  [band_num, data_num] = size(V);
  V = [V; delta*ones(1, data_num)]; % extern V for sum-to-one constraint
  
  % init W
  W = abs(randn(N, data_num));
  sum_W = sum(W);
  for i = 1:data_num
    W( :,i ) = W( :,i ) ./ sum_W(i);
  end
  
  % init H
  if HasInit == 1 % init H from nfidr
    H = V(1:band_num, H_InitIndx); % H(band_num, N)
  else
    H = abs( randn(band_num, N) );
  end
  H = [H; delta*ones(1, N)]; % H(N, N)
  H_record = zeros(band_num, N, IterMax);
  H_record(:, :, 1) = H(1:band_num, :);
  
  % calculate initial error
  e = (V - (H * W)) ./ V;
  e_sum = sum( sum( e.^2 ) );
  E = zeros(1, IterMax);
  
  % start iteration
  iter_cur = 1;
  E(iter_cur) = e_sum;
  mu = mean( V, 2); % get mean value for every band
  tau = lamda / factorial(N-1);
  alpha_k = alpha0;
  beta_k = beta0;
  while( e_sum > Tolerance)
    if(iter_cur == IterMax)
      break;
    end
    
    % the following update algorithm is according to the paper MVC
    % the symbol is the same for easy-understanding
    % update H
    U =  pca(V', 'Numcomponents', N-1); % (band_num, N-1)
    C = [ones(1, N); zeros(N-1, N)]; % (N, N)
    B = [zeros(1, N-1); eye(N-1, N-1)]; % (N, N-1)
    Z = C + B*U'*(H - mu*ones(1, N));
    Z_ = inv(Z);
    grad_H = (H*W - V)*W' + tau * det(Z)^2 * U * B' * Z_';
    H_prev = H;
    H = H - alpha_k * grad_H;
    ZEROS = zeros(band_num+1, N);
    H = max(H, ZEROS); % clip H great than 0
    
    % update step size for H
    f1 = mvc_object(V, H_prev, W, lamda);
    f2 = mvc_object(V, H, W, lamda);
    Delta_f = grad_H' * (H - H_prev);
    Delta_f = norm(Delta_f);
    mk = log( abs(f2-f1) / (sigma*alpha0*Delta_f) ) / log(tho);
    mk = floor(mk);
    mk = 0;
    alpha_k = tho^mk * alpha0;
    
    % update W
    grad_W = H' * (H*W - V);
    W_prev = W;
    W = W - beta_k * grad_W;
    ZEROS = zeros(N, data_num);
    W = max(W, ZEROS); % clip W great than 0
    
    % update step size for W
    f1 = f2;
    f2 = mvc_object(V, H, W, lamda);
    Delta_f = grad_W' * (W - W_prev);
    Delta_f = norm(Delta_f);
    mk = log( abs(f2-f1) / (sigma*beta0*Delta_f) ) / log(tho);
    mk = floor(mk);
    mk = 0;
    beta_k = tho^mk * beta0;
    
    % calculate SAD
    e = (V - (H * W)) ./ V;
    e_sum = sum( sum( e.^2 ) );
    
    iter_cur = iter_cur + 1;
    E(iter_cur) = e_sum;
    H_record(:, :, iter_cur) = H(1:band_num, :);
  end
  H_record = H_record(:, : , 1:iter_cur);
  E = E(1:iter_cur);
  % H = U' * (H - mu * ones(1, N));
  
end