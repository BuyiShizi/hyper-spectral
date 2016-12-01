function sum1 = nmfw(X,n)

sum1 = 0;

for i=1:n
    
    sum1 = sum1 + sum(X(i,:));
end
end