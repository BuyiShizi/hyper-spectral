%% this program create user-defined hyperspectrum data 
%  to verify the vadibility of the program for endmember extraction
clear all
%% create endmember spectrum data
%  we create 10 endmembers which cover 50 bands interval
endmember_num = 6;
band_num = 50;
endmember_data = zeros (endmember_num, band_num); %% endmember hyperspectrum data

band_data = 1:band_num;

endmember_data(1,:) = 4 * band_data;
endmember_data(2,:) = 150 * sin(band_data / band_num * 2) + 160;
endmember_data(3,:) = 150 * cos(band_data / band_num * 2) + 160;
endmember_data(4,:) = 150 * sin(band_data / band_num * 4) + 160;
endmember_data(5,:) = 150 * cos(band_data / band_num * 4) + 160;
endmember_data(6,:) = 150 * sin(band_data / band_num * 8) + 160;
figure(1)
plot (endmember_data');

%% create hyperspectrum data
data_num = 10000;
noise_level = 0;
hyper_data = ones (data_num, band_num);
abundance_fraction = zeros (data_num, endmember_num);

% for i = 1:10000
%     abundance_fraction(i,1:endmember_num-1) = rand(1,endmember_num-1);
%     abundance_fraction(i,endmember_num) = 1-sum(abundance_fraction(i,1:endmember_num-1));
% end
abundance_sub = zeros (10,endmember_num);
abundance_sub(1,:) = [1,   0,   0,   0,   0,   0];
abundance_sub(2,:) = [0,   1,   0,   0,   0,   0];
abundance_sub(3,:) = [0,   0,   1,   0,   0,   0];
abundance_sub(4,:) = [0,   0,   0,   1,   0,   0];
abundance_sub(5,:) = [0,   0,   0,   0,   1,   0];
abundance_sub(6,:) = [0,   0,   0,   0,   0,   1];
abundance_sub(7,:) = [0.1, 0.1, 0.15, 0.1, 0.45, 0.1];
abundance_sub(8,:) = [0.6, 0.1, 0.1, 0.1, 0.05 0.05];
abundance_sub(9,:) = [0.1, 0.5, 0.2, 0.2, 0, 0];
abundance_sub(10,:) = [0.1, 0.1, 0.1, 0.2, 0.4, 0.1];
abundance_fraction = repmat (abundance_sub, data_num/10, 1);
hyper_data = abundance_fraction * endmember_data;
hyper_data = hyper_data + noise_level*rand (data_num, band_num);
hyper_data_image = zeros(band_num, 100, 100);

for i = 1:band_num
    hyper_data_image(i,:,:) = reshape (hyper_data(:,i), 100, 100);
end

imshow_data = hyper_data_image(50,:,:);
imshow_data = reshape (imshow_data, 100, 100);
figure (2)
imshow (imshow_data/300);

%% endmember extraction
endmember_index = zeros (1,endmember_num);
[endmember_index, V] = n_findr (hyper_data, endmember_num);

%% plot the sprectrum of endmember
figure(3)
plot (hyper_data(endmember_index,:)')

%% calculate abundance
abundance_solution = zeros (endmember_num,1);
hyper_data_to_solve = hyper_data(1,:);
% abundance_solution = linear_least_square(hyper_data(endmember_index,:)', hyper_data_to_solve');
Y = hyper_data_to_solve;
% C = hyper_data(endmember_index,:);
C = endmember_data(1:6,:);
% X_INIT = abundance_solution';
% X_INIT = [0, 0, 0, 0, 0, 0];
X_INIT = rand(1,6);
X_INIT = X_INIT/sum(X_INIT)
X_INIT = [0 0 0 0 0 0];
RATE = 1e-3 ;
TOLERANT = 1e-6;
TOLERANT1 = 1;
%% plot error figure
% [x1,x2] = meshgrid (0:0.01:1,0:0.01:1);
% [x1_size, x2_size] = size(x1);
% e = zeros (x1_size, x2_size);
% for i = 1:x1_size
%     for j = 1:x2_size
%         f_error = Y' - C' * [x1(i,j); x2(i,j)];
%         e(i,j) = sum (f_error.^2);
%     end
% end
% contour(x1, x2, e,80); 


[abundance_solution1, final_error] = linear_least_square_gradient (Y', C', X_INIT', RATE, TOLERANT, TOLERANT1)
% abundance_solution = lsqnonneg(endmember_data', hyper_data(1,:)');