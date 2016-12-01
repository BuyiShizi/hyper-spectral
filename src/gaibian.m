load urban_162
%load matlab_true_4
load matlab_refer_nmf_one
Data=reshape(x,162,307,307);
 
 duanyuan=zeros(162,4);
 duanyuan(:,1)=Data(:,62,197);
 duanyuan(:,2)=Data(:,124,213);
 duanyuan(:,3)=Data(:,134,87);
 duanyuan(:,4)=Data(:,274,31);
O=zeros(307,307,4);
for i=1:307
    for j=1:307
        abundance= lsqnonneg(duanyuan,Data(:,i,j));%非负矩阵分解
        for k=1:4
            %if(abundance(k)>=0)
            O(i,j,k) = abundance(k);
          %  else if(abundance<=0)
              %      O(i,j,k)=0;
             %   end
           % end
        end
    end
end
%归一化
nnfengdu=reshape(O,94249,4);
 total=sum(nnfengdu,2);
  for j=1:94249
     nnfengdu(j,:)=nnfengdu(j,:)./total(j,1); 
  end
  fengdu_1=reshape(nnfengdu,307,307,4);
  %计算RMSE
  fengdu_1=double(fengdu_1);
                 S_est=double(S_est);
                 %s_est1=s_est(:,1);
                 MSE = sum(sum(abs(fengdu_1 - S_est).^2 ) );
                 [w,h]=size(fengdu_1);
                 TotalPixel = w*h;
                 RMSE=sqrt(MSE /TotalPixel);  
  %计算SAD
     
               duanyuan=double(duanyuan);
               A_est=double(A_est);
              % A_est1=A_est(:,1);
               %member=sum(a.*a_est,2);
               member=sum(duanyuan.*A_est);
               denominator=(sqrt(sum(A_est.^2,2)))'*(sqrt(sum(duanyuan.^2,2)));
               %sad= acos(member./denominator)
                [r,c]=size(member);
                %denominator=repmat(denominator,r,c);
               num=member./denominator;
               acos(num);