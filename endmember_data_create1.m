%% this program create user-defined hyperspectrum data 
%  to verify the vadibility of the program for endmember extraction
clear all
%% create endmember spectrum data
%  we create 10 endmembers which cover 50 bands interval
endmember_num = 6;
band_num = 50;
endmember_data = zeros (band_num, endmember_num); %% endmember hyperspectrum data

band_data = 1:band_num;

endmember_data(:,1) = 4 * band_data;
endmember_data(:,2) = 300 * sin(band_data / band_num * 2);
endmember_data(:,3) = 300 * cos(band_data / band_num * 2);
endmember_data(:,4) = 300 * sin(band_data / band_num * 4);
endmember_data(:,5) = 300 * cos(band_data / band_num * 4);
endmember_data(:,6) = 300 * sin(band_data / band_num * 8);
figure(1)
plot (endmember_data);

%% create hyperspectrum data
data_num = 10000;
noise_level = 0;
hyper_data = ones (band_num, data_num);
abundance_fraction = zeros (endmember_num, data_num);

% for i = 1:10000
%     abundance_fraction(i,1:endmember_num-1) = rand(1,endmember_num-1);
%     abundance_fraction(i,endmember_num) = 1-sum(abundance_fraction(i,1:endmember_num-1));
% end
abundance_sub = zeros (endmember_num,10);
abundance_sub(:,1) = [1,   0,   0,   0,   0,   0]';
abundance_sub(:,2) = [0,   1,   0,   0,   0,   0]';
abundance_sub(:,3) = [0,   0,   1,   0,   0,   0]';
abundance_sub(:,4) = [0,   0,   0,   1,   0,   0]';
abundance_sub(:,5) = [0,   0,   0,   0,   1,   0]';
abundance_sub(:,6) = [0,   0,   0,   0,   0,   1]';
abundance_sub(:,7) = [0.1, 0.1, 0.15, 0.1, 0.45, 0.1]';
abundance_sub(:,8) = [0.6, 0.1, 0.1, 0.1, 0.05 0.05]';
abundance_sub(:,9) = [0.1, 0.5, 0.2, 0.2, 0, 0]';
abundance_sub(:,10) = [0.1, 0.1, 0.1, 0.2, 0.4, 0.1]';
abundance_fraction = repmat (abundance_sub, 1,data_num/10);
hyper_data = endmember_data * abundance_fraction;
hyper_data = hyper_data + noise_level*rand (band_num, data_num);
hyper_data_image = zeros(band_num, 100, 100);

for i = 1:band_num
    hyper_data_image(i,:,:) = reshape (hyper_data(i,:), 100, 100);
end

imshow_data = hyper_data_image(50,:,:);
imshow_data = reshape (imshow_data, 100, 100);
figure (2)
imshow (imshow_data/300);

%% endmember extraction
endmember_index = zeros (1,endmember_num);
[endmember_index, V] = n_findr (hyper_data', endmember_num);

%% plot the sprectrum of endmember
figure(3)
plot (hyper_data(:,endmember_index))

%% calculate abundance
endmember_data_cal = hyper_data(:,endmember_index);
abundance_solution = zeros (endmember_num,data_num);

for i = 1:data_num
    abundance_solution(:,i) = linear_least_square(endmember_data_cal, hyper_data(:,i));
end

imshow_data = zeros(endmember_num, 100, 100);
for i = 1:endmember_num
    im=reshape(abundance_solution(i,:),100,100);
    imshow_data(i,:,:) = im;
end
figure (4)
for i = 1:endmember_num
    subplot (3,2,i);
    im = reshape(imshow_data(i,:,:),100,100);
    imshow (im);
end