%% create hyperdata 
%% the hyperdata consists two endmember data
clear all
%% create abundance map for two endmember data
clear all
figure_indx = 1;
hyper_data_size = 50;
x = linspace(-pi/2, pi/2, hyper_data_size);
y = linspace(-pi/2, pi/2, hyper_data_size);
[X, Y] = meshgrid(x, y);
level = 0;
Z1 = 1 * randn(hyper_data_size);
Z2 = 1 * randn(hyper_data_size);
Z3 = 1 * randn(hyper_data_size);
Z1 = abs(Z1);
Z2 = abs(Z2);
Z3 = abs(Z3);
Z_sum = Z1 + Z2 + Z3;
Z1 =  Z1 ./ Z_sum;
Z2 = Z2 ./ Z_sum;
Z3 = Z3 ./ Z_sum;
% load('Z1.mat')
% load('Z2.mat')
% load('Z3.mat')

% figure(1)
% subplot(131)
% mesh(X, Y, Z1); 
% subplot(132)
% mesh(X, Y, Z2); 
% subplot(133)
% mesh(X, Y, Z3); 
% mesh(X, Y, Z2);
% mesh(X, Y, Z3);

% figure(1)
% imshow(Z1); hold on
% imshow(Z2 + Z1)

%% create endmember data spectrum
spectrum_num = 10;
endmember_num  =  3;

f = 1:10;
e1 = rand(1,spectrum_num);
e2 = rand(1,spectrum_num);
e3 = rand(1,spectrum_num);
% load e1
% load e2
% load e3
% figure(2)
% plot(f, e1, f, e2, f, e3)


%% create hyperdata using endmember data spectrum and abundance map
hyper_data = zeros(spectrum_num, hyper_data_size, hyper_data_size);
hyper_data_1 = zeros(hyper_data_size, hyper_data_size);
hyper_data_sum = zeros(hyper_data_size, hyper_data_size);
hyper_data_3 = zeros(hyper_data_size, hyper_data_size,3);
noise_level = 0.01;
for i = 1:spectrum_num
    hyper_data_1 = e1(i) * Z1 + e2(i) * Z2 + e3(i) * Z3 + noise_level * rand(hyper_data_size);
    hyper_data_sum = hyper_data_1;
    hyper_data(i,:,:) = reshape(hyper_data_1, 1, hyper_data_size, hyper_data_size);
end

% figure(3)
% hyper_data_3(:,:,1) = reshape(hyper_data(1,:,:), hyper_data_size, hyper_data_size,1);
% hyper_data_3(:,:,2) = reshape(hyper_data(3,:,:), hyper_data_size, hyper_data_size,1);
% hyper_data_3(:,:,3) = reshape(hyper_data(5,:,:), hyper_data_size, hyper_data_size,1);
% imshow(hyper_data_3)

%% find endmember data using n_findr
alpha = 0.4;

hyper_data_line = zeros (hyper_data_size * hyper_data_size, spectrum_num);
for i = 1:spectrum_num
    hyper_data_line(:, i) = reshape(hyper_data(i, :, :), hyper_data_size * hyper_data_size, 1);
end
endmember_index = zeros (1,endmember_num);
[end_indx, V] = n_findr(hyper_data_line(:,1:2),3);
H_I = hyper_data_line( end_indx, 1:2 );
[W,H,H1, H2, H3, E,IT] = hyper_nmf ( hyper_data_line(:,1:2), 3, H_I, alpha, 0.001, 2000);
figure(figure_indx)
% subplot(221)
scatter( hyper_data_line(:,1), hyper_data_line(:,2) );
title('scatter view');
hold on
plot( hyper_data_line( end_indx(1),1 ), hyper_data_line( end_indx(1),2 ), 'ro')
plot( hyper_data_line( end_indx(2),1 ), hyper_data_line( end_indx(2),2 ), 'ro')
plot( hyper_data_line( end_indx(3),1 ), hyper_data_line( end_indx(3),2 ), 'ro')
plot( H1(1,:)/alpha, H1(2,:)/alpha,'c',H2(1,:)/alpha, H2(2,:)/alpha,'g',H3(1,:)/alpha, H3(2,:)/alpha,'y', 'LineWidth',2 )

% subplot(222)
% scatter( hyper_data_line(:,1), hyper_data_line(:,3) );
% title('scatter view');
% hold on
% plot( hyper_data_line( end_indx(1),1 ), hyper_data_line( end_indx(1),3 ), 'ro')
% plot( hyper_data_line( end_indx(2),1 ), hyper_data_line( end_indx(2),3 ), 'ro')
% plot( hyper_data_line( end_indx(3),1 ), hyper_data_line( end_indx(3),3 ), 'ro')
% plot( H1(1,:)/alpha, H1(3,:)/alpha,'c',H2(1,:)/alpha, H2(3,:)/alpha,'g',H3(1,:)/alpha, H3(3,:)/alpha,'y', 'LineWidth',2 )
% 
% subplot(223)
% scatter( hyper_data_line(:,2), hyper_data_line(:,3) );
% title('scatter view');
% hold on
% plot( hyper_data_line( end_indx(1),2 ), hyper_data_line( end_indx(1),3 ), 'ro')
% plot( hyper_data_line( end_indx(2),2 ), hyper_data_line( end_indx(2),3 ), 'ro')
% plot( hyper_data_line( end_indx(3),2 ), hyper_data_line( end_indx(3),3 ), 'ro')
% plot( H1(2,:)/alpha, H1(3,:)/alpha,'c',H2(2,:)/alpha, H2(3,:)/alpha,'g',H3(2,:)/alpha, H3(3,:)/alpha,'y', 'LineWidth',2 )
%% test
test_line = W * H;
figure(5)
hold on
h_scatter = scatter(test_line(:,endm_ind1)/alpha, test_line(:,endm_ind2)/alpha, 'r');
%% rotate hyper data 
theta = 270;
T = [cosd(theta), -sind(theta), 0; sind(theta), cos(theta), 0 ; 0, 0, 1];
hyper_data_line_T = hyper_data_line(:,1:3) * T;
scatter(hyper_data_line_T(:,1), hyper_data_line_T(:,2));
hold on
H_T = H(:,1:3) * T;
plot ( H_T(1,endm_ind1)/alpha, H_T(1,endm_ind2)/alpha, 'ko');
plot ( H_T(2,endm_ind1)/alpha, H_T(2,endm_ind2)/alpha, 'ko');
plot ( H_T(3,endm_ind1)/alpha, H_T(3,endm_ind2)/alpha, 'ko');
% H123_T = [H1',H2',H3'] * inv(T);
H1_T = H1' * T(1:2,1:2);
H2_T = H2' * T(1:2,1:2);
H3_T = H3' * T(1:2,1:2);
plot(H1_T(:,1)/alpha,H1_T(:,2)/alpha, 'k');
plot(H2_T(:,1)/alpha,H2_T(:,2)/alpha, 'k');
plot(H3_T(:,1)/alpha,H3_T(:,2)/alpha, 'k');
%%
[endmember_index, V, M] = n_findr_orig (hyper_data_line, endmember_num);
figure(4)
plot (hyper_data_line(endmember_index,:)')
M1 = M(1,:);
M2 = M(2,:);
M3 = M(3,:);
M_index1 = find(M(1,:)>0);
M_index2 = find(M(2,:)>0);
M_index3 = find(M(3,:)>0);

%%
scatter_index = 1:hyper_data_size * hyper_data_size;
endm_ind1 = 1;
endm_ind2 = 2;
record_start = 1;
record_end = IT-1;
figure(5)
hold on
h_scatter = scatter(hyper_data_line(:,endm_ind1), hyper_data_line(:,endm_ind2));
% plot ( H(1,endm_ind1)/alpha, H(1,endm_ind2)/alpha, 'ko');
% plot ( H(2,endm_ind1)/alpha, H(2,endm_ind2)/alpha, 'ko');
% plot ( H(3,endm_ind1)/alpha, H(3,endm_ind2)/alpha, 'ko');
% plot(H1(1,record_start:record_end)/alpha,H1(2,record_start:record_end)/alpha, 'k');
% plot(H2(1,record_start:record_end)/alpha,H2(2,record_start:record_end)/alpha, 'k');
% plot(H3(1,record_start:record_end)/alpha,H3(2,record_start:record_end)/alpha, 'k');
% index1_size = length(M_index1);
% annotation('arrow',[0.1, 0.2],[0.1, 0.2], 'Color', 'r')
% m1_x = hyper_data_line(M1(M_index1),endm_ind1);
% m1_y = hyper_data_line(M1(M_index1),endm_ind2);
% plot_arrow(m1_x,m1_y, 'r', 2)
% 
% m2_x = hyper_data_line(M2(M_index2),endm_ind1);
% m2_y = hyper_data_line(M2(M_index2),endm_ind2);
% plot_arrow(m2_x,m2_y, 'g', 2)
% 
% m3_x = hyper_data_line(M3(M_index3),endm_ind1);
% m3_y = hyper_data_line(M3(M_index3),endm_ind2);
% plot_arrow(m3_x,m3_y, 'k', 2)
%% solve abundance map using normal linear square method
E = [e1'; e2']; % endmember data
abundance_map_line = zeros (hyper_data_size * hyper_data_size, endmember_num);

for i = 1:hyper_data_size * hyper_data_size
    hyper_data_to_solve = hyper_data_line(i, :);
    abundance_solution = linear_least_square(hyper_data_line(endmember_index,:)', hyper_data_to_solve');
    abundance_map_line(i, :) = abundance_solution;
end
 
abundance_map_1 = reshape(abundance_map_line(:,1), hyper_data_size, hyper_data_size);
abundance_map_2 = reshape(abundance_map_line(:,2), hyper_data_size, hyper_data_size);
abundance_map_3 = reshape(abundance_map_line(:,3), hyper_data_size, hyper_data_size);

figure(6)
hold on
subplot(311)
imshow(abundance_map_1)
subplot(312)
imshow(abundance_map_2)
subplot(313)
imshow(abundance_map_3)
% imshow(abundance_map_2)
% imshow(abundance_map_1)

%% solve abundance map using gradient method
% E = [e1', e2'];
% abundance_map = zeros (hyper_data_size * hyper_data_size, 2);
% 
% X_INIT = zeros (hyper_data_size * hyper_data_size, 2);
% RATE = 1e-3 ;
% TOLERANT = 1e-6;
% [X, ERROR] = gradient_descent(hyper_data_line', E, X_INIT', RATE, TOLERANT);
% 
% abundance_map_grad_1 = reshape (X(1,:), hyper_data_size, hyper_data_size);
% abundance_map_grad_2 = reshape (X(2,:), hyper_data_size, hyper_data_size);
% figure(7)
% imshow(abundance_map_grad_2);
abundance_map_nmf_1 = 0;
abundance_map_nmf_2 = 0;
abundance_map_nmf_3 = 0;
%% solve abundance map nonnegative solution woth NMF-similiar method
E = [e1', e2', e3'];
Tolerance = 0.05;
MaxIter = 10000;
[X, error] = gradient_descent_nonnegative(hyper_data_line', E, MaxIter, Tolerance);

% abundance_map_nmf_1_before = abundance_map_nmf_1;
% abundance_map_nmf_2_before = abundance_map_nmf_2;
% abundance_map_nmf_3_before = abundance_map_nmf_3;
abundance_map_nmf_1 = reshape (X(1,:), hyper_data_size, hyper_data_size);
abundance_map_nmf_2 = reshape (X(2,:), hyper_data_size, hyper_data_size);
abundance_map_nmf_3 = reshape (X(3,:), hyper_data_size, hyper_data_size);
figure(8)
hold on
subplot(231)
imshow(abundance_map_nmf_1);
subplot(232)
imshow(abundance_map_nmf_2);
subplot(233)
imshow(abundance_map_nmf_3);
% subplot(234)
% imshow(abundance_map_nmf_1_before)
% subplot(235)
% imshow(abundance_map_nmf_2_before)
% subplot(236)
% imshow(abundance_map_nmf_3_before)
% imshow(abundance_map_nmf_2);
% imshow(abundance_map_nmf_1);
% 
% %% create hyperdata using endmember data spectrum and abundance map
% hyper_data_nmf = zeros(spectrum_num, hyper_data_size, hyper_data_size);
% hyper_data_1 = zeros(hyper_data_size, hyper_data_size);
% hyper_data_sum = zeros(hyper_data_size, hyper_data_size);
% hyper_data_3 = zeros(hyper_data_size, hyper_data_size,3);
% noise_level = 0;
% for i = 1:spectrum_num
%     hyper_data_1 = e1(i) * abundance_map_nmf_1 + e2(i) * abundance_map_nmf_2 + e3(i) * abundance_map_nmf_3 + noise_level * rand(hyper_data_size);
%     hyper_data_sum = hyper_data_1;
%     hyper_data_nmf(i,:,:) = reshape(hyper_data_1, 1, hyper_data_size, hyper_data_size);
% end
% 
% hyper_data_3(:,:,1) = reshape(hyper_data_nmf(1,:,:), hyper_data_size, hyper_data_size,1);
% hyper_data_3(:,:,2) = reshape(hyper_data_nmf(3,:,:), hyper_data_size, hyper_data_size,1);
% hyper_data_3(:,:,3) = reshape(hyper_data_nmf(5,:,:), hyper_data_size, hyper_data_size,1);
% 
% hyper_data_line_nmf = zeros (hyper_data_size * hyper_data_size, spectrum_num);
% for i = 1:spectrum_num
%     hyper_data_line_nmf(:, i) = reshape(hyper_data_nmf(i, :, :), hyper_data_size * hyper_data_size, 1);
% end
% [endmember_index, V] = n_findr_orig (hyper_data_line_nmf, endmember_num);
% 
% % before
% hyper_data_nmf_before = zeros(spectrum_num, hyper_data_size, hyper_data_size);
% hyper_data_1 = zeros(hyper_data_size, hyper_data_size);
% hyper_data_sum = zeros(hyper_data_size, hyper_data_size);
% hyper_data_3 = zeros(hyper_data_size, hyper_data_size,3);
% noise_level = 0;
% for i = 1:spectrum_num
%     hyper_data_1 = e1(i) * abundance_map_nmf_1_before + e2(i) * abundance_map_nmf_2_before + e3(i) * abundance_map_nmf_3_before + noise_level * rand(hyper_data_size);
%     hyper_data_sum = hyper_data_1;
%     hyper_data_nmf_before(i,:,:) = reshape(hyper_data_1, 1, hyper_data_size, hyper_data_size);
% end
% 
% hyper_data_3(:,:,1) = reshape(hyper_data_nmf_before(1,:,:), hyper_data_size, hyper_data_size,1);
% hyper_data_3(:,:,2) = reshape(hyper_data_nmf_before(3,:,:), hyper_data_size, hyper_data_size,1);
% hyper_data_3(:,:,3) = reshape(hyper_data_nmf_before(5,:,:), hyper_data_size, hyper_data_size,1);
% 
% hyper_data_line_nmf_before = zeros (hyper_data_size * hyper_data_size, spectrum_num);
% for i = 1:spectrum_num
%     hyper_data_line_nmf_before(:, i) = reshape(hyper_data_nmf_before(i, :, :), hyper_data_size * hyper_data_size, 1);
% end
% [endmember_index_before, V] = n_findr_orig (hyper_data_line_nmf_before, endmember_num);
% 
% % 
% figure(9)
% hold on
% scatter(hyper_data_line_nmf(:,endm_ind1), hyper_data_line_nmf(:,endm_ind2))
% scatter(hyper_data_line_nmf(endmember_index,endm_ind1),hyper_data_line_nmf(endmember_index,endm_ind2),'r','filled');
% scatter(hyper_data_line_nmf_before(endmember_index_before,endm_ind1),hyper_data_line_nmf_before(endmember_index_before,endm_ind2),'k','filled');