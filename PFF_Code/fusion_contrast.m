function y = fusion_contrast(tv,tr,xv,xr,yv,yr,a,ir,ir_b)
% tv,tr--待融合的可见光与红外对比度图像
% xv,xr--结构显著性信息
% yv,yr--像素级信息
% a--分解层数
% ir--红外图像
%ir_b--红外多尺度分解图像

%% 红外权重图调整系数

[hei, wid] = size(ir);

sr = ir-ir_b;

p = 1./(1+exp((0.03).*(0-sr)));

for i = 1 : hei
        for j = 1 : wid
            if p(i,j) >=0.5
                p(i,j)=2*p(i,j) ;
            else
                p(i,j)=2*p(i,j)/(a);
            end 
        end
end

%% 逐像素显著性和结构显著性
P_V = abs(yv);  P_R = abs(yr); % 逐像素显著性

r=[2,6,12,24,48,96];

S_V = gradient_f(xv,r(a)); S_R = gradient_f(xr,r(a)); %结构显著性

wv = 1*P_V+5*S_V;  wr = p.*(1*P_R+5*S_R); %显著性合并

d = wv-wr;

W = 1./(1+exp((0.1*2.^(4-a)).*(0-d))); %计算权重

c = a; 
b1 = fspecial('gaussian',[3*c+1,3*c+1],c/2);

W = imfilter(W,b1,'replicate');

%% 融合

y = W.*tv+(1-W).*tr;  %加权融合


end