function F = transducer(R,t)
% R-����ֲ���Ӧ�ԶԱȶ�
% F-�����֪�Աȶ�
% t-�м����

w = 22;
u =0.5; %��ֵ
[m,n] = size(R);
F = ones(m,n);
J = ones(m,n);
J(R<=0) = -1;
R = abs(R);
for i = 1:m
    for j = 1:n
            if R(i,j)>= u
                F(i,j) = w*(R(i,j)/u)^(1/2);
            else
                F(i,j) = w*(R(i,j)/u)^t;
            end
    end
end
F=F.*J;
end





