%% Verion 1
% This progrom reshape scatter plot first,
%  and then, do nmf-endmember parsing
hyper_data_shape_prev = hyper_data_line(:,1:3);
hyper_data_shape_cur    = hyper_data_shape_prev(:,1:3)';
AT_A = hyper_data_shape_cur' * hyper_data_shape_cur;
load('hyper_data_ref.mat');
hyper_data_trans = hyper_data_ref;
clear hyper_data_ref
tmp = hyper_data_shape_cur' * hyper_data_trans;
trans_matrix = AT_A \ tmp;
hyper_data_shape_aft = hyper_data_shape_cur * trans_matrix;

hyper_data_shape_nmf = hyper_data_shape_aft';
[end_indx_shape, V_shape, M_shape] = n_findr_orig(hyper_data_shape_nmf,3);
H_I = hyper_data_shape_nmf( end_indx_shape, 1:3 );
[W_shape,H_shape,H1_shape, H2_shape, H3_shape, E_shape,IT_shape] = hyper_nmf ( hyper_data_shape_nmf, 3, H_I, alpha, 0.001, 2000);
figure_index = 1;
figure(figure_index)
figure_index = figure_index + 1;
subplot(221)
scatter( hyper_data_shape_nmf(:,1), hyper_data_shape_nmf(:,2) );
title('scatter view');
hold on
plot( hyper_data_shape_nmf( end_indx_shape(1),1 ), hyper_data_shape_nmf( end_indx_shape(1),2 ), 'ro')
plot( hyper_data_shape_nmf( end_indx_shape(2),1 ), hyper_data_shape_nmf( end_indx_shape(2),2 ), 'ro')
plot( hyper_data_shape_nmf( end_indx_shape(3),1 ), hyper_data_shape_nmf( end_indx_shape(3),2 ), 'ro')
plot( H1_shape(1,:)/alpha, H1_shape(2,:)/alpha,'c',H2_shape(1,:)/alpha, H2_shape(2,:)/alpha,'g',H3_shape(1,:)/alpha, H3_shape(2,:)/alpha,'y', 'LineWidth',2 )

subplot(222)
scatter( hyper_data_shape_nmf(:,1), hyper_data_shape_nmf(:,3) );
title('scatter view');
hold on
plot( hyper_data_shape_nmf( end_indx_shape(1),1 ), hyper_data_shape_nmf( end_indx_shape(1),3 ), 'ro')
plot( hyper_data_shape_nmf( end_indx_shape(2),1 ), hyper_data_shape_nmf( end_indx_shape(2),3 ), 'ro')
plot( hyper_data_shape_nmf( end_indx_shape(3),1 ), hyper_data_shape_nmf( end_indx_shape(3),3 ), 'ro')
plot( H1_shape(1,:)/alpha, H1_shape(3,:)/alpha,'c',H2_shape(1,:)/alpha, H2_shape(3,:)/alpha,'g',H3_shape(1,:)/alpha, H3_shape(3,:)/alpha,'y', 'LineWidth',2 )

subplot(223)
scatter( hyper_data_shape_nmf(:,2), hyper_data_shape_nmf(:,3) );
title('scatter view');
hold on
plot( hyper_data_shape_nmf( end_indx_shape(1),2 ), hyper_data_shape_nmf( end_indx_shape(1),3 ), 'ro')
plot( hyper_data_shape_nmf( end_indx_shape(2),2 ), hyper_data_shape_nmf( end_indx_shape(2),3 ), 'ro')
plot( hyper_data_shape_nmf( end_indx_shape(3),2 ), hyper_data_shape_nmf( end_indx_shape(3),3 ), 'ro')
plot( H1_shape(2,:)/alpha, H1_shape(3,:)/alpha,'c',H2_shape(2,:)/alpha, H2_shape(3,:)/alpha,'g',H3_shape(2,:)/alpha, H3_shape(3,:)/alpha,'y', 'LineWidth',2 )

 recorver original scatter
trans_matrix_inv = inv(trans_matrix);
h1 = H1_shape(:,2000)/alpha;
h2 = H2_shape(:,2000)/alpha;
h3 = H3_shape(:,2000)/alpha;
hyper_data_shape_nmf( end_indx_shape,:) = [h1';h2';h3'];
hyper_data_shape_nmf_rec = hyper_data_shape_nmf' * trans_matrix_inv;

figure(figure_index)
figure_index = figure_index + 1;
subplot(121)
scatter( hyper_data_shape_nmf(:,1), hyper_data_shape_nmf(:,2))
hold on
plot( hyper_data_shape_nmf( end_indx_shape(1),1 ), hyper_data_shape_nmf( end_indx_shape(1),2 ), 'ro')
plot( hyper_data_shape_nmf( end_indx_shape(2),1 ), hyper_data_shape_nmf( end_indx_shape(2),2 ), 'ro')
plot( hyper_data_shape_nmf( end_indx_shape(3),1 ), hyper_data_shape_nmf( end_indx_shape(3),2 ), 'ro')
axis([0 1 0 1])
axis equal
subplot(122)
scatter( hyper_data_shape_nmf_rec(1,:), hyper_data_shape_nmf_rec(2,:))
hold on
plot( hyper_data_shape_nmf_rec( 1,end_indx_shape(1) ), hyper_data_shape_nmf_rec( 2,end_indx_shape(1) ), 'ro')
plot( hyper_data_shape_nmf_rec( 1,end_indx_shape(2) ), hyper_data_shape_nmf_rec( 2,end_indx_shape(2) ), 'ro')
plot( hyper_data_shape_nmf_rec( 1,end_indx_shape(3) ), hyper_data_shape_nmf_rec( 2,end_indx_shape(3) ), 'ro')
axis([0 1 0 1])
axis equal