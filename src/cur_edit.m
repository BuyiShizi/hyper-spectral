%% simple triangle trasfrom
a = [ 0.3 0.3; 0.5 0.35; 1 0.35]';
a_mean = mean(a');
a_mean = repmat( a_mean', 1, 3);
a = a - a_mean;

C = a * a';
[U, S, V] = svd( C );
M = inv(sqrt(S)) * V';
b = M * a;

plot( a(1,:), a(2,:), 'g')
hold on
plot( [a(1,1), a(1,3)],[a(2,1),a(2,3)],'g')


plot( b(1,:), b(2,:), 'r')
plot( [b(1,1), b(1,3)],[b(2,1),b(2,3)],'r')
axis equal

%% random triangle transfrom
a = rand( 2, 50 )
C = a * a';
[U, S, V] = svd( C );
M = inv(sqrt(S)) * V';
b = M * a;

scatter( a(1,:), a(2,:) )
hold on
scatter( b(1,:), b(2,:) , 'r')
%% real triangle transfrom
figure_index = 1;
point_set = hyper_data_line( :,1:2 );
point_set = point_set';
AT_A = point_set' * point_set;
% AT_A_inv = inv( AT_A );
load('hyper_data_ref.mat');
point_set_trans = hyper_data_ref;
clear hyper_data_ref
tmp  = point_set' * point_set_trans;
trans_matrix = AT_A \ tmp;
point_set_trans_s = point_set *  trans_matrix;

figure(figure_index)
figure_index = figure_index + 1;
subplot( 131 )
scatter( point_set_trans(1,:), point_set_trans(2,:) );
axis([0 1 0 1])
axis equal
hold on
scatter( point_set_trans_s(1,:), point_set_trans_s(2,:),  'r' );

subplot( 132 )
scatter( point_set_trans_s(1,:), point_set_trans_s(2,:))
axis([0 1 0 1])
axis equal

subplot( 133 )
scatter( point_set(1,:), point_set(2,:))
axis([0 1 0 1])
axis equal

% recorver original scatter
trans_matrix_inv = inv(trans_matrix);
point_set_rec = point_set_trans_s * trans_matrix_inv; 
figure(figure_index)
figure_index = figure_index + 1;
subplot(121)
scatter( point_set(1,:), point_set(2,:))
axis([0 1 0 1])
axis equal
subplot(122)
scatter( point_set_rec(1,:), point_set_rec(2,:))
axis([0 1 0 1])
axis equal

