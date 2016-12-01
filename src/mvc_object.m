function f = mvc_object(V, H, W, lamda)
%% object function in nmf_MVC algorithm
% V        : original hyper data matrix[band_num, data_size]
% H,W      : current factorized matrix,
%            H[band_num, unit_num]
%            W[unit_num, data_size]
% lamda    : factor to control the tradeoff between
%            the accurate reconstruction and the volumn constraint

%% start process
  % calculate least linear square
  [band_num, unit_num] = size(H);
  lls = norm(V - H*W, 'fro')^2 / 2;
 
  % calculate the cost function
  H_ = H - mean(V,2)*ones(1, unit_num);
  pc = pca(V', 'Numcomponents', unit_num-1); % pca analyis
  H_ = pc' * H_;
  H_ = [ones(1, unit_num); H_];
  J = det(H_)^2 / (2*factorial(unit_num-1));
  
  f = lls + lamda * J;