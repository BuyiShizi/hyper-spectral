function N_FINDR()

% read a hyperspectral image from sample data of ENVI.
nlines  = 307; %254
npixels = 307; %3454
nbands  = 210;
nlength = nlines*npixels; % 每个波段读出来的像素的个数
%hs_img  = multibandread('cup95eff.int', [nlines npixels nbands], ...
                       % 'int16', 0, 'bil', 'ieee-le'); %从二进制文件中读取BSQ,BIL,BIP数据
%hs_img= multibandread('F:\285-destripe_new_flaash_wv_sg_37',[nlines npixels nbands],...
 %                     'double',0,'bsq','ieee-be');%从二进制文件中读取BSQ文件
% test_hdr='F:\test_img.hdr';%头文件的地址
% maindata = ReadImg(test_hdr);%读取头文件相应的.img高光谱数据

S=fopen('URBAN2','r');
% M=307;
% N=307;
% band=210;
R=zeros(nlines,npixels,nbands);
for i=1:nbands
    F=fread(S,[nlines npixels],'int16');
    R(:,:,i)=F';
end
hs_img=R;

% reduct the dimensionality of image to be one less the number of end-
% mumbers by using PCA.
% standardize the data by dividing each column by its standard deviation.
hs_img = reshape(hs_img, [], nbands);%改变矩阵的形状，但是元素个数不变 data_num_per_band * nbands
meanh  = mean(hs_img);%求均值
stdh   = std(hs_img);%求标准差
sd_img = (hs_img-repmat(meanh,nlength,1))./repmat(stdh,nlength,1);%repmat：将meanh作为元素复制nlength*1块%标准化图像 data_num_per_band * nbands

%做PCA。~是对主分的打分也就是原矩阵在主成分空间的表示，pcoef是元矩阵对应的协方差矩阵对的所有特征向量，latent协方差矩阵的特征值
[pcoef, score, latent] = princomp(sd_img);
                                       

% determine the number of endmembers by calculating the contribution of
% %计算主成分的贡献来确定端元的数量
% principal components.
perc = cumsum(latent)/sum(latent)*100;%cumsum计算一个数据各行的累加值返回值与latent行列相同，sum是将latent的所有元素相加
Nend = sum(perc<99.5)+1;%确定端元的数量
% get the principal components.确定主成分的贡献量
pca_img = sd_img*pcoef(:,1:Nend-1); % data_num_per_band * main_nbands

% repeat 50 times to find the endmembers with largest
% volume.%迭代50次找到最大体积的端元
Ntimes = 50;
locs   = zeros(Ntimes,Nend);
V_max  = zeros(Ntimes,1);
for i = 1:Ntimes%tic和toc计算程序的时间
   
    [locs(i,:), V_max(i,1)] = finder(pca_img, nlength, Nend);%N-finder的函数
   
end
[score, ind] = max(V_max);
loc = locs(ind, :);

% calculate the abundance of each endmember for each
% pixel.计算对于每个像素的每个端元的丰度（最小二乘法）
M   = ones(Nend);
M(2:end,:) = pca_img(loc,:)';
C = zeros(Nend, nlength);
for i = 1:nlength
    p = [1; pca_img(i,:)'];
    C(:,i) = lsqnonneg(M, p);%返回C>=0,约束下norm(M*C-p)的最小值向量
%     C(:,i) = lsqlin(E, p, [],[],[],[],0,1);
end

save('F:\NFINDR_res.txt','C','-ASCII');

%stop;

% end function N_FINDR.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [loc, V_max] = finder(pca_img, nlength, Nend)

% randomly select pixels as initial endmembers.随机选择像素作为初始端元，求单形体的体积
ind = unidrnd(nlength, 1, Nend);%产生从1到nlength所指定的最大数数之间的离散均匀随机数，生成1*nend矩阵 1 * Nend
E   = ones(Nend);%产生nend*nend全是1的数组 Nend * Nend
E(2:end,:) = pca_img(ind,:)';%随机选择nend个像素组成最初的端元矩阵
dentor = factorial(Nend-1);%求（端元-1）的阶乘
V_max  = abs(det(E))/dentor;%求E的行列式在取绝对值最后除以dentor。求体积

% find the largest volume, and set the corresponding pixels as
% endmembers.找到最大的体积，并设置相应的像素作为端元
for i = 1:Nend  % loop over each endmember.循环每个端元
    i_max = ind(i);
    for j = 1:nlength  % loop over each pixel.循环每个像素
        E(2:end,i) = pca_img(j,:)';
        V = abs(det(E))/dentor;
        k = V>V_max;
        i_max(k) = j;
        V_max(k) = V;
    end
    E(2:end,i) = pca_img(i_max,:)';
    ind(i) = i_max;
end
loc = ind;

% end function finder().
end


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 读取img格式文件
% 分两步进行：1.解析头文件 2.按照文件信息读取数据
% 返回数据矩阵
% 解析头文件
% 输入hdr文件
% 输出信息：图像宽，高，层数，数据类型，排列类型，光谱信息

function maindata = ReadImg(hdrfilename)
filename=hdrfilename;
fidhdr=fopen(filename);
hdrdata=textscan(fidhdr,'%s');
fclose(fidhdr);

hdrmatrix=hdrdata{1,1};
[m,n]=size(hdrmatrix);
%  samples = num
  samplesNumIndex=FindStrIndex(hdrmatrix,'samples')+2;
%   str2num转换时要进行强制转换
        hdrInfo.samples=str2num(char(hdrmatrix(samplesNumIndex)));
%  lines = num
          linesNumIndex=FindStrIndex(hdrmatrix,'lines')+2;
        hdrInfo.lines=str2num(char(hdrmatrix(linesNumIndex)));
%  bands = num
         bandsNumIndex= FindStrIndex(hdrmatrix,'bands')+2;
         hdrInfo.bands=str2num(char(hdrmatrix(bandsNumIndex)));
%         默认offset为零  data type = num
datatypeNumIndex= FindStrIndex(hdrmatrix,'data')+3;
hdrInfo.datatype=str2num(char(hdrmatrix(datatypeNumIndex)));
%  interleave = '***'
interleaveIndex=FindStrIndex(hdrmatrix,'interleave')+2;
hdrInfo.interleave= char(hdrmatrix(interleaveIndex));
disp(hdrInfo.interleave);


% 根据头文件信息读取数据
% 主要区别在于数据类型和数据排列方式,只写出5种常用数据类型
% 数据类型
 datatypestr={'uint8','int16','int32','single','double'};
% datatypestr={'uint8','int32','single','double'};
% 排布方式,统一 lower到小写比较
interleaveStyles={'bsq','bil','bip'};
 
fidimg=fopen(strrep(filename,'hdr','img')); 
% 清除data
clear data;
switch hdrInfo.interleave
%     隔帧存储
    case interleaveStyles(1)
        for i=1: hdrInfo.bands
            data(:,:,i)=fread(fidimg,[hdrInfo.samples,hdrInfo.lines],char(datatypestr(hdrInfo.datatype)));
        end    
%      隔行存储
    case interleaveStyles(2)
        for i=1:hdrInfo.lines
            for j=1:hdrInfo.bands
                  data(i,:,j)=fread(fidimg,hdrInfo.samples,char(datatypestr(hdrInfo.datatype)));
            end
        end
%         隔点存储
    case  interleaveStyles(3)
         for i=1:hdrInfo.lines
                 for j=1:hdrInfo.samples
                     for k=1:hdrInfo.bands
                      data(i,j,k)=fread(fidimg,1,char(datatypestr(hdrInfo.datatype)));
                     end
                 end
         end
end
fclose(fidimg); 
maindata=data;
 WriteImg('test.img',maindata,hdrInfo.datatype);
end
%  写img 文件
% 分别写头文件和数据文件
% 默认bsq方式写入
function WriteImg(filename,dataMatrix,datatype)
datatypestr={'uint8','int16','int32','single','double'};
% 由矩阵得出文件信息
sizeInfo=size(dataMatrix);
hdr.lines=sizeInfo(1);
hdr.samples=sizeInfo(2);
if length(sizeInfo)==3
hdr.bands=sizeInfo(3);
else
    hdr.bands=1;
end
hdr.datatype=datatypestr(datatype);
% 准备文件名
if isempty(strfind(filename,'.img'))
   hdrname= strcat(filename,'.hdr');
else
   hdrname=strrep(filename,'.img','.hdr')
end
% 一项一项的输入
fidhdr=fopen(hdrname,'w');
fprintf(fidhdr,'%s\r\n','ENVI'); 
fprintf(fidhdr,'description = {\r\nCreated Time: %s}\r\n',date);
fprintf(fidhdr,'samples = %d\r\n',hdr.samples);
fprintf(fidhdr,'lines = %d\r\n',hdr.lines);
fprintf(fidhdr,'bands = %d\r\n',hdr.bands);
fprintf(fidhdr,'%s\r\n','file type = ENVI Standard');
fprintf(fidhdr,'data type = %d\r\n',datatype);
 fprintf(fidhdr,'%s\r\n','interleave = bsq');
 fprintf(fidhdr,'%s\r\n','sensor type = Unknown');
 fprintf(fidhdr,'%s\r\n','byte order = 0');
 fprintf(fidhdr,'%s\r\n','wavelength units = Unknown');
%  暂时不加波段描述信息
%  fprintf(fidhdr,'%s\r\n','band names = Unknown'); 
fclose(fidhdr);

% 写入数据文件
% 准备文件名
if isempty(strfind(filename,'.img'))
   imgname= strcat(filename,'.img'); 
else
    imgname=filename;
end
imgfid=fopen(imgname,'w');
for k=1:hdr.bands 
%     for i=1:hdr.lines
%         for j=1:hdr.samples
%          fwrite(imgfid,dataMatrix,hdr.datatype);
%         end
%     end
    fwrite(imgfid,dataMatrix(:,:,k),char(hdr.datatype));
end  
fclose(imgfid);
end


% 查找字符串的位置
function  Index=FindStrIndex(hdrMatrix,objstr)
[m,n]=size(hdrMatrix);
for i=1:m
    if strcmp(hdrMatrix(i),objstr)==1        
        break;
    end
end
Index=i;
end



