%% 读入源图像
%=================================================================================================
% close all;
global rpos;
rpos = rotatePosition(80); 

is_RGB = false;   % The setting when the visible image is a gray image
% is_RGB = true;  % The setting when the visible image is a RGB color image
 z = 1; % Select the No.of source image pair in ".\images\Gray" or ".\images\RGB" depending on is_RGB 

if is_RGB
    % RGB image pairs stored in ".\images\RGB"
    pv=['.\images\RGB\VIS\',num2str(z),'.jpg'];           
    pr=['.\images\RGB\IR\',num2str(z),'.jpg'];
    [img_V,Cb,Cr] =  image_imread_rgb_v(pv);
    img_R = image_imread_rgb_r(pr);
    figure;imshow(imread(pv)); title('Visible image');
    figure;imshow(img_R); title('Infrared image');
else
    % Gray image pairs stored in ".\images\Gray"
    pv=['.\images\Gray\VIS\',num2str(z),'.png'];           
    pr=['.\images\Gray\IR\',num2str(z),'.png']; 
    img_V = image_imread(pv);
    img_R = image_imread(pr);
    figure;imshow(img_V); title('Visible image');
    figure;imshow(img_R); title('Infrared image');
end
disp(['Conducting the fusion ...'])
%=================================================================================================
%% 定义变量及参数
%=================================================================================================
nLevel =4; %分解层数
t = [1.40,1.15,1.04,1.15,1.35,1.93]; 
I_V = cell(1,nLevel+1);    I_R = cell(1,nLevel+1); %多尺度分解
M_V = cell(1,nLevel);      M_R = cell(1,nLevel); %局部对比度函数
D_V = cell(1,nLevel);      D_R = cell(1,nLevel); %带通图像
R_V =  cell(1,nLevel);     R_R =  cell(1,nLevel);%局部带通对比度图像
V_V = cell(1,nLevel);      V_R = cell(1,nLevel); %非线性转换对比度图像
C_V = cell(1,nLevel);      C_R = cell(1,nLevel); %去噪以及抑制过饱和后对比度图像
C_F = cell(1,nLevel);  %融合对比度图像
%=================================================================================================

%% 多尺度分解
%=================================================================================
%非加速版本
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I_V{1} = img_V;  I_R{1} = img_R;
% lambda1 =0.1;
% I_V{2} = SAEP(I_V{1}, lambda1,4);I_R{2} = SAEP(I_R{1}, lambda1,4);
% lambda2 =1;
% I_V{3} = SAEP(I_V{2}, lambda2,8,4);I_R{3} = SAEP(I_R{2}, lambda2,8,4);
% lambda3 =1.9;
% I_V{4} = SAEP(I_V{3}, lambda3,16,8);I_R{4} = SAEP(I_R{3}, lambda3,16,8);
% lambda4 =2.8;
% I_V{5} = SAEP(I_V{4}, lambda4,32,16);I_R{5} = SAEP(I_R{4}, lambda4,32,16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 加速版本
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_V{1} = img_V;  I_R{1} = img_R;
lambda1 =0.1;
I_V{2} = SAEP(I_V{1}, lambda1,4);I_R{2} = SAEP(I_R{1}, lambda1,4);
lambda2 =1;
I_V{3} = SAEP(I_V{2}, lambda2,8,4);I_R{3} = SAEP(I_R{2}, lambda2,8,4);
I_V{3}= imresize(I_V{3},0.5,'bicubic');I_R{3}= imresize(I_R{3},0.5,'bicubic');
lambda3 =1.9;
I_V{4} = SAEP(I_V{3}, lambda3,8,4);I_R{4} = SAEP(I_R{3}, lambda3,8,4);
lambda4 =2.8;
I_V{5} = SAEP(I_V{4}, lambda4,16,8);I_R{5} = SAEP(I_R{4}, lambda4,16,8);

I_V{3}= imresize(I_V{3},2,'bicubic');I_R{3}= imresize(I_R{3},2,'bicubic');
I_V{4}= imresize(I_V{4},2,'bicubic');I_R{4}= imresize(I_R{4},2,'bicubic');
I_V{5}= imresize(I_V{5},2,'bicubic');I_R{5}= imresize(I_R{5},2,'bicubic');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=================================================================================
%% 人类视觉响应空间转换及噪声与强度饱和抑制
%===========================================================================
P_V = 0; P_R = 0; %聚合初始化
for i = 1 : nLevel
    D_V{i} =255.* (I_V{i} -I_V{i+1}); D_R{i} =255.*(I_R{i} - I_R{i+1}); 

    M_V{i} = gainv(255.*I_V{i+1}); M_R{i} = gainr(255.*I_R{i+1}); 

    R_V{i} = M_V{i}.*D_V{i}; R_R{i} =M_R{i}.*D_R{i}; 

    V_V{i} = transducer(R_V{i},t(i)); V_R{i} = transducer(R_R{i},t(i));

    %噪声和强度饱和抑制
    C_V{i} = detajf(V_V{i},I_V{1},i); C_R{i} = detajf(V_R{i},I_R{1},i);

    P_V = P_V+C_V{i}; P_R = P_R+C_R{i};
end
% =========================================================================
%% 基础层图像处理及融合
%==========================================================================
B_V = sbgain(255.*I_V{nLevel+1});   B_R = sbgain(255.*I_R{nLevel+1});  

T_V = B_V+C_V{nLevel}; T_R = B_R+C_R{nLevel};

c=16; 
s_v = saliency(T_V,c*3); s_r = saliency(T_R,c*3);
w_v = s_v./(s_v+s_r);
b = fspecial('gaussian',[6*c+1,6*c+1],c);
w_v = imfilter(w_v,b,'replicate');

B_F = w_v.*B_V+(1-w_v).*B_R;
%==========================================================================
%% 带通对比度图像融合
%==========================================================================
S_V = B_V+C_V{nLevel};S_R = B_R+C_R{nLevel};

for i = nLevel:-1:2
    if i<nLevel
        P_V = P_V-C_V{i+1}; P_R = P_R-C_R{i+1};  %自上而下像素级信息
        S_V = S_V+C_V{i};  S_R = S_R+C_R{i}; %自下而上结构显著性
    end

    C_F{i} = fusion_contrast(C_V{i},C_R{i},S_V,S_R,P_V,P_R,i,img_R.*255,I_R{i+1}.*255);
end

C_F{1} = fusion_contrast(C_V{1},C_R{1},S_V,S_R,C_V{1},C_R{1},1,img_R.*255,I_R{2}.*255);

%========================================================================================================
%% 图像重构
%================================================
for i = nLevel:-1:1  
    %带通对比度inv-transducer
    R_F = in_transducer(C_F{i},t(i));

    %带通
    AL = R_F./gainv(B_F);
  
    %相加得到上一层次低通图像
    B_F = AL + B_F;

end
%================================================
%% 图像显示
%==========================================================================
if is_RGB
    img_y = B_F; 
    YCbCr1(:,:,1)=uint8(img_y);
    YCbCr1(:,:,2)=Cb;
    YCbCr1(:,:,3)=Cr;
    img_F = ycbcr2rgb(YCbCr1);
    figure;imshow(img_F); 
%     imwrite(img_F,['.\results\RGB\',num2str(z),'.1','.png']);
    clear YCbCr1
else
    img_F = B_F; 
    img0 = uint8(img_F);
    figure;imshow(img0); title('Fused image');  
%     imwrite(img_F/255,['.\results\Gray\',num2str(z),'.1','.png']);
end
disp(['Finish!'])
% ========================================================================