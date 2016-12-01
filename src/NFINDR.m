tic
clear all;
clc
im=multibandread('D:\1result\mnf',[200 200 5],'float32',0,'bil','ieee-le',...
{'Band','Range',[1 4]});   %读取高光谱图像数据
A=reshape(im,200*200,4);%重构光谱向量
A=A';%A矩阵的每一列表示一个像素的光谱向量
p=randperm(200*200);%随机生成数
k=p(1:5);%取p中前5个数
B(:,:)=zeros(4,5);%初始化B矩阵
results=zeros(1,5);%初始化results向量
for i=1:5   %随机选择5个像元作为端元的初始向量
B(:,i)=A(:,k(i));
end
E0=[1,1,1,1,1;B];
V0=abs(det(E0))/24;   %计算初始端元的单形体体积
Vmax=V0;
for i=1:5           %逐步替换端元向量，找出使单形体体积最大的端元向量
    for j=1:40000
        B(:,i)=A(:,j);
        Eij=[1,1,1,1,1;B];
        Vij=abs(det(Eij))/24;
        if Vij>Vmax            
            Emax=[1,1,1,1,1;B];
            Vmax=Vij; 
            results(i)=j;            
        end
    end
    B(:,i)=A(:,results(i));
end
location(:,:)=zeros(5,2);
for i=1:5
    temp=num2location(results(i),200,200);
    location(i,:)=temp;
end
disp(location);
toc

