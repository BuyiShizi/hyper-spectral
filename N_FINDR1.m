function N_FINDR()

% read a hyperspectral image from sample data of ENVI.
nlines  = 307; %254
npixels = 307; %3454
nbands  = 210;
nlength = nlines*npixels; % ÿ�����ζ����������صĸ���
%hs_img  = multibandread('cup95eff.int', [nlines npixels nbands], ...
                       % 'int16', 0, 'bil', 'ieee-le'); %�Ӷ������ļ��ж�ȡBSQ,BIL,BIP����
%hs_img= multibandread('F:\285-destripe_new_flaash_wv_sg_37',[nlines npixels nbands],...
 %                     'double',0,'bsq','ieee-be');%�Ӷ������ļ��ж�ȡBSQ�ļ�
% test_hdr='F:\test_img.hdr';%ͷ�ļ��ĵ�ַ
% maindata = ReadImg(test_hdr);%��ȡͷ�ļ���Ӧ��.img�߹�������

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
hs_img = reshape(hs_img, [], nbands);%�ı�������״������Ԫ�ظ������� data_num_per_band * nbands
meanh  = mean(hs_img);%���ֵ
stdh   = std(hs_img);%���׼��
sd_img = (hs_img-repmat(meanh,nlength,1))./repmat(stdh,nlength,1);%repmat����meanh��ΪԪ�ظ���nlength*1��%��׼��ͼ�� data_num_per_band * nbands

%��PCA��~�Ƕ����ֵĴ��Ҳ����ԭ���������ɷֿռ�ı�ʾ��pcoef��Ԫ�����Ӧ��Э�������Ե���������������latentЭ������������ֵ
[pcoef, score, latent] = princomp(sd_img);
                                       

% determine the number of endmembers by calculating the contribution of
% %�������ɷֵĹ�����ȷ����Ԫ������
% principal components.
perc = cumsum(latent)/sum(latent)*100;%cumsum����һ�����ݸ��е��ۼ�ֵ����ֵ��latent������ͬ��sum�ǽ�latent������Ԫ�����
Nend = sum(perc<99.5)+1;%ȷ����Ԫ������
% get the principal components.ȷ�����ɷֵĹ�����
pca_img = sd_img*pcoef(:,1:Nend-1); % data_num_per_band * main_nbands

% repeat 50 times to find the endmembers with largest
% volume.%����50���ҵ��������Ķ�Ԫ
Ntimes = 50;
locs   = zeros(Ntimes,Nend);
V_max  = zeros(Ntimes,1);
for i = 1:Ntimes%tic��toc��������ʱ��
   
    [locs(i,:), V_max(i,1)] = finder(pca_img, nlength, Nend);%N-finder�ĺ���
   
end
[score, ind] = max(V_max);
loc = locs(ind, :);

% calculate the abundance of each endmember for each
% pixel.�������ÿ�����ص�ÿ����Ԫ�ķ�ȣ���С���˷���
M   = ones(Nend);
M(2:end,:) = pca_img(loc,:)';
C = zeros(Nend, nlength);
for i = 1:nlength
    p = [1; pca_img(i,:)'];
    C(:,i) = lsqnonneg(M, p);%����C>=0,Լ����norm(M*C-p)����Сֵ����
%     C(:,i) = lsqlin(E, p, [],[],[],[],0,1);
end

save('F:\NFINDR_res.txt','C','-ASCII');

%stop;

% end function N_FINDR.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [loc, V_max] = finder(pca_img, nlength, Nend)

% randomly select pixels as initial endmembers.���ѡ��������Ϊ��ʼ��Ԫ������������
ind = unidrnd(nlength, 1, Nend);%������1��nlength��ָ�����������֮�����ɢ���������������1*nend���� 1 * Nend
E   = ones(Nend);%����nend*nendȫ��1������ Nend * Nend
E(2:end,:) = pca_img(ind,:)';%���ѡ��nend�������������Ķ�Ԫ����
dentor = factorial(Nend-1);%�󣨶�Ԫ-1���Ľ׳�
V_max  = abs(det(E))/dentor;%��E������ʽ��ȡ����ֵ������dentor�������

% find the largest volume, and set the corresponding pixels as
% endmembers.�ҵ������������������Ӧ��������Ϊ��Ԫ
for i = 1:Nend  % loop over each endmember.ѭ��ÿ����Ԫ
    i_max = ind(i);
    for j = 1:nlength  % loop over each pixel.ѭ��ÿ������
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
% ��ȡimg��ʽ�ļ�
% ���������У�1.����ͷ�ļ� 2.�����ļ���Ϣ��ȡ����
% �������ݾ���
% ����ͷ�ļ�
% ����hdr�ļ�
% �����Ϣ��ͼ����ߣ��������������ͣ��������ͣ�������Ϣ

function maindata = ReadImg(hdrfilename)
filename=hdrfilename;
fidhdr=fopen(filename);
hdrdata=textscan(fidhdr,'%s');
fclose(fidhdr);

hdrmatrix=hdrdata{1,1};
[m,n]=size(hdrmatrix);
%  samples = num
  samplesNumIndex=FindStrIndex(hdrmatrix,'samples')+2;
%   str2numת��ʱҪ����ǿ��ת��
        hdrInfo.samples=str2num(char(hdrmatrix(samplesNumIndex)));
%  lines = num
          linesNumIndex=FindStrIndex(hdrmatrix,'lines')+2;
        hdrInfo.lines=str2num(char(hdrmatrix(linesNumIndex)));
%  bands = num
         bandsNumIndex= FindStrIndex(hdrmatrix,'bands')+2;
         hdrInfo.bands=str2num(char(hdrmatrix(bandsNumIndex)));
%         Ĭ��offsetΪ��  data type = num
datatypeNumIndex= FindStrIndex(hdrmatrix,'data')+3;
hdrInfo.datatype=str2num(char(hdrmatrix(datatypeNumIndex)));
%  interleave = '***'
interleaveIndex=FindStrIndex(hdrmatrix,'interleave')+2;
hdrInfo.interleave= char(hdrmatrix(interleaveIndex));
disp(hdrInfo.interleave);


% ����ͷ�ļ���Ϣ��ȡ����
% ��Ҫ���������������ͺ��������з�ʽ,ֻд��5�ֳ�����������
% ��������
 datatypestr={'uint8','int16','int32','single','double'};
% datatypestr={'uint8','int32','single','double'};
% �Ų���ʽ,ͳһ lower��Сд�Ƚ�
interleaveStyles={'bsq','bil','bip'};
 
fidimg=fopen(strrep(filename,'hdr','img')); 
% ���data
clear data;
switch hdrInfo.interleave
%     ��֡�洢
    case interleaveStyles(1)
        for i=1: hdrInfo.bands
            data(:,:,i)=fread(fidimg,[hdrInfo.samples,hdrInfo.lines],char(datatypestr(hdrInfo.datatype)));
        end    
%      ���д洢
    case interleaveStyles(2)
        for i=1:hdrInfo.lines
            for j=1:hdrInfo.bands
                  data(i,:,j)=fread(fidimg,hdrInfo.samples,char(datatypestr(hdrInfo.datatype)));
            end
        end
%         ����洢
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
%  дimg �ļ�
% �ֱ�дͷ�ļ��������ļ�
% Ĭ��bsq��ʽд��
function WriteImg(filename,dataMatrix,datatype)
datatypestr={'uint8','int16','int32','single','double'};
% �ɾ���ó��ļ���Ϣ
sizeInfo=size(dataMatrix);
hdr.lines=sizeInfo(1);
hdr.samples=sizeInfo(2);
if length(sizeInfo)==3
hdr.bands=sizeInfo(3);
else
    hdr.bands=1;
end
hdr.datatype=datatypestr(datatype);
% ׼���ļ���
if isempty(strfind(filename,'.img'))
   hdrname= strcat(filename,'.hdr');
else
   hdrname=strrep(filename,'.img','.hdr')
end
% һ��һ�������
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
%  ��ʱ���Ӳ���������Ϣ
%  fprintf(fidhdr,'%s\r\n','band names = Unknown'); 
fclose(fidhdr);

% д�������ļ�
% ׼���ļ���
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


% �����ַ�����λ��
function  Index=FindStrIndex(hdrMatrix,objstr)
[m,n]=size(hdrMatrix);
for i=1:m
    if strcmp(hdrMatrix(i),objstr)==1        
        break;
    end
end
Index=i;
end



