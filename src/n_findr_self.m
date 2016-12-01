%% read a hyperspectral image from sample data of ENVI.
clc;
clear all;

nlines  = 307; %254
npixels = 307; %3454
nbands  = 210;
nlength = nlines*npixels; % number of pixel in one band
%hs_img  = multibandread('cup95eff.int', [nlines npixels nbands], ...
                       % 'int16', 0, 'bil', 'ieee-le'); 
%hs_img= multibandread('F:\285-destripe_new_flaash_wv_sg_37',[nlines npixels nbands],...
 %                     'double',0,'bsq','ieee-be');
% test_hdr='F:\test_img.hdr';
% maindata = ReadImg(test_hdr);

S=fopen('URBAN2','r');
% M=307;
% N=307;
% band=210;
R=zeros(nlines,npixels,nbands);
for i=1:nbands
    F=fread(S,[nlines npixels],'int16');
    R(:,:,i) = F';
end
hs_img = R;

%% reduct the dimensionality of image to be one less the number of end-
% members by using PCA.
% standardize the data by dividing each column by its standard deviation.
hs_img = reshape(hs_img, [], nbands);% transform 2-dementional data into one dementional vector data_num_per_band * nbands
meanh  = mean(hs_img);% calculate average value for every band
stdh   = std(hs_img);% calculate std value for every band
sd_img = (hs_img-repmat(meanh,nlength,1))./repmat(stdh,nlength,1); % data_num_per_band * nbands

% PCA analysis
[pcoef, score, latent] = princomp(sd_img);
                                       

% determine the number of endmembers by calculating the contribution of
% 
% principal components.
perc = cumsum(latent)/sum(latent)*100;%
Nend = sum(perc<99.5)+1;
% get the principal components.
pca_img = sd_img*pcoef(:,1:Nend-1); % data_num_per_band * main_nbands

%% repeat 50 times to find the endmembers with largest
% volume.%迭代50次找到最大体积的端元
Ntimes = 50;
locs   = zeros(Ntimes,Nend);
V_max  = zeros(Ntimes,1);
for i = 1:Ntimes % number of iteration to do N_FINDT algorithm
   disp(i);
    [locs(i,:), V_max(i,1)] = finder(pca_img, nlength, Nend); %N-finder function
   
end
[score, ind] = max(V_max); % get the most satisfying result
loc = locs(ind, :);

%% plot the curve for endmember spectrum
figure
hold on
title ('endmember spectrun');
xlabel ('band');
ylabel ('value');

band_index = 1:nbands;
endmember_data = zeros (1,nbands);
color = ['r','g','k','m'];
for i=1:4
    endmember_data = hs_img(loc(i),:);
    plot (band_index, endmember_data,color(i));
end

%% calculate the abundance of each endmember for each endmember
% the method for canculate abundance is linear least square with nonnegative constraints with fully constraints:
% 1:the sum of pecentage abundance is 1
% 2:abundance coefficient is great than 0
M   = ones(Nend);
M(2:end,:) = pca_img(loc,:)';
C = zeros(Nend, nlength);
for i = 1:nlength
    disp(i);
    p = [1; pca_img(i,:)'];
    C(:,i) = lsqnonneg(M, p); %穕inear least square with nonnegative constraints, it is a full constraints
%     C(:,i) = lsqlin(E, p, [],[],[],[],0,1);
end


%% plot the abundance of the hyspectral
map_n = 4; % show how many bands for abundance map
abundance_fractional = C(1,:);
abundance_map = zeros (map_n,nlength);
for (i=1:map_n)
	abundance_map(i,:) = abundance_fractional .* hs_img(:,i)';
end
abundance_map_reshape = reshape (round(abundance_map(1,:)), 307, 307);
figure (2)
title('abundance map');
imshow(abundance_map_reshape);