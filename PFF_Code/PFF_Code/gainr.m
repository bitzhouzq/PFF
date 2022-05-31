function F =  gainr(M)
Mmin = min(M(:));
if Mmin < 0 
    M = M-Mmin;
end
c=0.8; 
m=mean(M(:))/2;

F = 1./(1*(M+1).^c)+0.00;
[a,b]=size(F);

for i=1:a
    for j=1:b
        if M(i,j)<m
            F(i,j) = 1./(1*(M(i,j)+30).^c)+0.00;
        end
    end
end
end