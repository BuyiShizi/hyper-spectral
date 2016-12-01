tic
clear all;
clc
im=multibandread('D:\1result\mnf',[200 200 5],'float32',0,'bil','ieee-le',...
{'Band','Range',[1 4]});   %��ȡ�߹���ͼ������
A=reshape(im,200*200,4);%�ع���������
A=A';%A�����ÿһ�б�ʾһ�����صĹ�������
p=randperm(200*200);%���������
k=p(1:5);%ȡp��ǰ5����
B(:,:)=zeros(4,5);%��ʼ��B����
results=zeros(1,5);%��ʼ��results����
for i=1:5   %���ѡ��5����Ԫ��Ϊ��Ԫ�ĳ�ʼ����
B(:,i)=A(:,k(i));
end
E0=[1,1,1,1,1;B];
V0=abs(det(E0))/24;   %�����ʼ��Ԫ�ĵ��������
Vmax=V0;
for i=1:5           %���滻��Ԫ�������ҳ�ʹ������������Ķ�Ԫ����
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

