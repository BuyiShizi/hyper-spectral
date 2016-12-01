clc;
close all;
clear all;

%im1=imread('C:\Users\Administrator\Desktop\remote sensing\Reno_Radiance_wIGMGLT\band1.tif');
im1 = multibandread('hyperspectral_data.dat', [145 145 220],  'uint16', 0, 'bsq', 'ieee-le');
for i=185:190
imb=im1(:,:,i);
figure;imshow(imb,[])
end
% imb=imb-min(min(imb));
% imb=imb/max(max(imb));
%imb=double(imb);
