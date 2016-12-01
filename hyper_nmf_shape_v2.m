%% reshape endmember using Whitening

% whiten 1st and 2st data
hyper_unit = hyper_data_line( end_indx, 1:2)';
mean_unit = mean( hyper_unit, 2 );
mean_unit = repmat( mean_unit, 1, 3);
hyper_unit = hyper_unit - mean_unit;

C = hyper_unit * hyper_unit';
[U, S, V] = svd( C );
M = inv(sqrt(S)) * V';
hyper_unit_shape = M * hyper_unit;
ind1 = 1;
ind2 = 2;
plot( hyper_unit(ind1, :), hyper_unit(ind2, :), 'g')
hold on
plot( [hyper_unit( ind1, 1 ), hyper_unit( ind1, 3 )],[hyper_unit( ind2, 1 ),hyper_unit( ind2, 3 )],'g')

plot( hyper_unit_shape( ind1,:), hyper_unit_shape( ind2,:), 'r')
plot( [hyper_unit_shape( ind1,1), hyper_unit_shape( ind1,3)],[hyper_unit_shape( ind2, 1 ),hyper_unit_shape( ind2, 3 )],'r')
axis equal


%% solver again
    hyper_data_shape = hyper_data_line(:,1:2)';
    mean_hyper = mean( hyper_data_shape, 2);
    % mean_hyper = repmat( mean_hyper, 1, hyper_data_size*hyper_data_size);
    hyper_data_shape = [ hyper_data_shape; ones(1, hyper_data_size*hyper_data_size)];
    M_mean = [ 1 0 -mean_hyper(1); 0 1 -mean_hyper(2); 0 0 1];
    hyper_data_shape = M_mean * hyper_data_shape;
    hyper_data_shape_M  = M * hyper_data_shape(1:2,:);
    mini_hyper = min( hyper_data_shape_M, [], 2 );
    mini_hyper = repmat( mini_hyper, 1, hyper_data_size*hyper_data_size );
    scale = 1.1;
    hyper_data_shape_M = [ hyper_data_shape_M; ones(1,hyper_data_size*hyper_data_size)];
    M_shft = [ 1 0 -scale*mini_hyper(1); 0 1 -scale*mini_hyper(2); 0 0 1];
    hyper_data_shape_M = M_shft * hyper_data_shape_M;
    % hyper_data_shape_M  = hyper_data_shape_M - scale*mini_hyper;

% solve new scatter
    [end_indx_shape, V_shape] = n_findr(hyper_data_shape_M(1:2,:)',3);
    H_I = hyper_data_shape_M( 1:2,end_indx_shape );
    IT_max = 2000;
    [W_shape,H_shape,H1_shape, H2_shape, H3_shape, E_shape,IT_shape] = hyper_nmf ( hyper_data_shape_M(1:2,:)', 3, H_I', alpha, 0.001, IT_max);

% recover original scatters
    H1_recon = [ H1_shape(1:2,:)/alpha; ones(1, IT_max)];
    H2_recon = [ H2_shape(1:2,:)/alpha; ones(1, IT_max)];
    H3_recon = [ H3_shape(1:2,:)/alpha; ones(1, IT_max)];
    H1_recon = inv(M_shft) * H1_recon;
    H2_recon = inv(M_shft) * H2_recon;
    H3_recon = inv(M_shft) * H3_recon;
    H1_recon = inv(M) * H1_recon(1:2,:);
    H2_recon = inv(M) * H2_recon(1:2,:);
    H3_recon = inv(M) * H3_recon(1:2,:);
    H1_recon = [ H1_recon; ones(1, IT_max)];
    H2_recon = [ H2_recon; ones(1, IT_max)];
    H3_recon = [ H3_recon; ones(1, IT_max)];
    H1_recon = inv(M_mean) * H1_recon;
    H2_recon = inv(M_mean) * H2_recon;
    H3_recon = inv(M_mean) * H3_recon;
    hyper_data_recon = [hyper_data_shape_M];
    hyper_data_recon = inv(M_shft) * hyper_data_recon;
    hyper_data_recon = inv(M) * hyper_data_recon(1:2,:);
    hyper_data_recon = [hyper_data_recon; ones(1, hyper_data_size * hyper_data_size)];
    hyper_data_recon = inv(M_mean) * hyper_data_recon;
    
 % visualize
    figure_index = 1;
    figure(figure_index)
    figure_index = figure_index + 1;
    subplot(131)
    scatter( hyper_data_shape_M(1,:), hyper_data_shape_M(2,:) );
    title('scatter view');
    hold on
    axis equal
    plot( hyper_data_shape_M( 1, end_indx_shape(1) ), hyper_data_shape_M( 2, end_indx_shape(1) ), 'ro')
    plot( hyper_data_shape_M( 1, end_indx_shape(2) ), hyper_data_shape_M( 2, end_indx_shape(2) ), 'ro')
    plot( hyper_data_shape_M( 1, end_indx_shape(3) ), hyper_data_shape_M( 2, end_indx_shape(3) ), 'ro')
    plot( H1_shape(1,:)/alpha, H1_shape(2,:)/alpha,'c',H2_shape(1,:)/alpha, H2_shape(2,:)/alpha,'g',H3_shape(1,:)/alpha, H3_shape(2,:)/alpha,'y', 'LineWidth',2 )

    subplot(132)
    scatter( hyper_data_line(:,1), hyper_data_line(:,2) );
    title('scatter view');
    hold on
    % axis equal
    plot( hyper_data_line( end_indx(1),1 ), hyper_data_line( end_indx(1),2 ), 'ro')
    plot( hyper_data_line( end_indx(2),1 ), hyper_data_line( end_indx(2),2 ), 'ro')
    plot( hyper_data_line( end_indx(3),1 ), hyper_data_line( end_indx(3),2 ), 'ro')
    plot( H1(1,:)/alpha, H1(2,:)/alpha,'c',H2(1,:)/alpha, H2(2,:)/alpha,'g',H3(1,:)/alpha, H3(2,:)/alpha,'y', 'LineWidth',2 );

    subplot(133)
    scatter( hyper_data_recon(1,:), hyper_data_recon(2,:) );
    title('scatter view');
    hold on
%     axis equal
    plot( hyper_data_recon( 1,end_indx_shape(1) ), hyper_data_recon( 2,end_indx_shape(1) ), 'ro')
    plot( hyper_data_recon( 1,end_indx_shape(2) ), hyper_data_recon( 2,end_indx_shape(2) ), 'ro')
    plot( hyper_data_recon( 1,end_indx_shape(3) ), hyper_data_recon( 2,end_indx_shape(3) ), 'ro')
    plot( H1_recon(1,:), H1_recon(2,:),'c',H2_recon(1,:), H2_recon(2,:),'g',H3_recon(1,:), H3_recon(2,:),'y', 'LineWidth',2 );
    
% plot( hyper_data_shape_nmf( end_indx_shape(1),1 ), hyper_data_shape_nmf( end_indx_shape(1),3 ), 'ro')
% plot( hyper_data_shape_nmf( end_indx_shape(2),1 ), hyper_data_shape_nmf( end_indx_shape(2),3 ), 'ro')
% plot( hyper_data_shape_nmf( end_indx_shape(3),1 ), hyper_data_shape_nmf( end_indx_shape(3),3 ), 'ro')
% plot( H1_shape(1,:)/alpha, H1_shape(3,:)/alpha,'c',H2_shape(1,:)/alpha, H2_shape(3,:)/alpha,'g',H3_shape(1,:)/alpha, H3_shape(3,:)/alpha,'y', 'LineWidth',2 )
% 
% subplot(223)
% scatter( hyper_data_shape_nmf(:,2), hyper_data_shape_nmf(:,3) );
% title('scatter view');
% hold on
% plot( hyper_data_shape_nmf( end_indx_shape(1),2 ), hyper_data_shape_nmf( end_indx_shape(1),3 ), 'ro')
% plot( hyper_data_shape_nmf( end_indx_shape(2),2 ), hyper_data_shape_nmf( end_indx_shape(2),3 ), 'ro')
% plot( hyper_data_shape_nmf( end_indx_shape(3),2 ), hyper_data_shape_nmf( end_indx_shape(3),3 ), 'ro')
% plot( H1_shape(2,:)/alpha, H1_shape(3,:)/alpha,'c',H2_shape(2,:)/alpha, H2_shape(3,:)/alpha,'g',H3_shape(2,:)/alpha, H3_shape(3,:)/alpha,'y', 'LineWidth',2 )
% 
%  recorver original scatter
% trans_matrix_inv = inv(trans_matrix);
% h1 = H1_shape(:,2000)/alpha;
% h2 = H2_shape(:,2000)/alpha;
% h3 = H3_shape(:,2000)/alpha;
% hyper_data_shape_nmf( end_indx_shape,:) = [h1';h2';h3'];
% hyper_data_shape_nmf_rec = hyper_data_shape_nmf' * trans_matrix_inv;
% 
% figure(figure_index)
% figure_index = figure_index + 1;
% subplot(121)
% scatter( hyper_data_shape_nmf(:,1), hyper_data_shape_nmf(:,2))
% hold on
% plot( hyper_data_shape_nmf( end_indx_shape(1),1 ), hyper_data_shape_nmf( end_indx_shape(1),2 ), 'ro')
% plot( hyper_data_shape_nmf( end_indx_shape(2),1 ), hyper_data_shape_nmf( end_indx_shape(2),2 ), 'ro')
% plot( hyper_data_shape_nmf( end_indx_shape(3),1 ), hyper_data_shape_nmf( end_indx_shape(3),2 ), 'ro')
% axis([0 1 0 1])
% axis equal
% subplot(122)
% scatter( hyper_data_shape_nmf_rec(1,:), hyper_data_shape_nmf_rec(2,:))
% hold on
% plot( hyper_data_shape_nmf_rec( 1,end_indx_shape(1) ), hyper_data_shape_nmf_rec( 2,end_indx_shape(1) ), 'ro')
% plot( hyper_data_shape_nmf_rec( 1,end_indx_shape(2) ), hyper_data_shape_nmf_rec( 2,end_indx_shape(2) ), 'ro')
% plot( hyper_data_shape_nmf_rec( 1,end_indx_shape(3) ), hyper_data_shape_nmf_rec( 2,end_indx_shape(3) ), 'ro')
% axis([0 1 0 1])
% axis equal