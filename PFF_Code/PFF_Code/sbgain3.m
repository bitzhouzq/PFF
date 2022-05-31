function y = sbgain3(x,a)
% x-低通输入图像
% y-低通适应性处理输出
[m,n] = size(x);
y = zeros(m,n);

for i=1:m
    for j =1:n
        if x(i,j)<128 %低于128
            if a<3 
                y(i,j) = x(i,j)*gainv(x(i,j))/gainv(128);
            else
                y(i,j) = x(i,j)*gainv(x(i,j))/gainv(x(i,j));
            end
        else
            if a<1
                y(i,j) = x(i,j)*gainv(x(i,j))/gainv(128);
            else 
                y(i,j) = x(i,j)*gainv(x(i,j))/gainv(x(i,j));
            end
        end
    end
end
t = 1;
b = fspecial('gaussian',[6*t+1,6*t+1],t);
% b = fspecial('gaussian',[7,7],1);
y = imfilter(y,b,'replicate');
end