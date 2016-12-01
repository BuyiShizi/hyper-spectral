S=fopen('URBAN2','r');
M=307;
N=307;
band=210;
R=zeros(M,N,band);
for i=1:band
    F=fread(S,[M N],'int16');
    R(:,:,i)=F';
end
RR=zeros(M,N,162);
for i=1:71
    RR(:,:,i)=R(:,:,i+4);
end
for i=72:81
    RR(:,:,i)=R(:,:,i+5);
end
for i=82:94
    RR(:,:,i)=R(:,:,i+6);
end
for i=95:118
    RR(:,:,i)=R(:,:,i+17);
end
for i=119:162
    RR(:,:,i)=R(:,:,i+35);
end
X = zeros(162,M,N);
TT = zeros(M,N);
for i=1:162
    % TT = RR(:,:,i);
    for j=1:M
        for l=1:N
            TT(j,l)=RR(j,l,i);
        end
    end
    X(i,:)=TT(:)';
    X(i,:)=X(i,:)/sum(X(i,:));
end
W=nmfw(X,4);
y=W*X;
figure 
imshow(reshape(y(1,:),307,307),[])
figure 
imshow(reshape(y(50,:),307,307),[])
figure 
imshow(reshape(y(100,:),307,307),[])
figure 
imshow(reshape(y(150,:),307,307),[])


