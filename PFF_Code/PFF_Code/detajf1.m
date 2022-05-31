function y = detajf1(x,g,l) %噪声和过饱和的抑制
% x--输入图像
% g--参考图像，取源图像
% l--第l层图像
[m,n] = size(x);
y1 = x;


% 抑制过饱和
for i = 1:m
    for j = 1:n
        if x(i,j)>0&&l>2
            y1(i,j) = x(i,j)*(g(i,j).^1.5*0.7+0.3);
%         elseif x(i,j)<0&&l>2                                 %*********Zhou**
%             y1(i,j) = x(i,j)*((1.01-g(i,j)).^1.5*0.7+0.3);   %*********Zhou**

        end
    end
end
y=y1;

for i = 1:m
    for j = 1:n
        if y1(i,j)>50
            y(i,j) =50;
        else
            if y1(i,j)< -50
            y(i,j) = -50;
            end

        end
    end
end


%抑制噪声

if l==1
    down = ((255-mean(255.*g(:)))/255)^2*2.^(-l+3); 
    y(abs(x)<down) = 0;

else
    if l<3 && l>1
    down = ((255-mean(255.*g(:)))/255)^2*2;
    y(abs(x)<down) = 0;
    end
end
end