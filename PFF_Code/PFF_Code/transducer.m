function F = transducer(M,a)
% M-输入局部适应性对比度
% F-输出感知对比度
% a-中间参数

t1 = 22;
t0 = t1;
t2 =0.5; %阈值，小于阈值抑制
[m,n] = size(M);
F = ones(m,n);
for i = 1:m
    for j = 1:n
        if M(i,j)>0 
            if M(i,j)>t2
                F(i,j) = t1*(M(i,j)/t2)^(1/2);
            else
                F(i,j) = t1*(M(i,j)/t2)^a;
            end
        else 
            M(i,j) = abs(M(i,j));
            if M(i,j)>t2
                F(i,j) = -1*t0*(M(i,j)/t2)^(1/2);
            else
                F(i,j) = -1*t0*(M(i,j)/t2)^a;
            end
        end
    end
end


end





