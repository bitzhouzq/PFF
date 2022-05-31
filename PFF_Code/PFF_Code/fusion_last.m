function y = fusion_last(tv,tr,xv,xr,yv,yr,a,ir,b,bb)
% tv,tr--待融合可见光与红外视觉对比度
% xv,xr--结构显著性叠加
% yv,yr--（细节显著性）
% ir--红外图像

%% 求红外图像的显著图

[hei, wid] = size(ir);

sr = ir-b;

s =1./(1+exp((0.03).*(0-sr))); %显著图

for i = 1 : hei
        for j = 1 : wid
            if s(i,j) >=0.5
                s(i,j)=2*s(i,j) ;
            else
                s(i,j)=2*s(i,j)/(a);
            end 
        end
end


%% 逐像素显著性和结构显著性
dsv = abs(yv);  dsr = abs(yr); % 逐像素显著性

r=[2,6,12,24];
ssv = gradient_f(xv,r(a)); ssr = gradient_f(xr,r(a)); %结构显著性

wv = 1*dsv+5*ssv;  wr = s.*(1*dsr+5*ssr); %显著性合并

d = wv-wr;

w = 1./(1+exp((0.1*2.^(4-a)).*(0-d))); %计算权重

c = a; 
b1 = fspecial('gaussian',[3*c+1,3*c+1],c/2);

w1 = imfilter(w,b1,'replicate');

%% 融合

y = w1.*tv+(1-w1).*tr;  %加权融合


end