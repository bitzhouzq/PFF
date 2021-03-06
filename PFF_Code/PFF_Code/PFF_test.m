%% 源图像
%================================================================================
close all;
global rpos;
rpos = rotatePosition(80);
 for z=1
   p1=['.\images\VIS\',num2str(z),'.png'];           
   p2=['.\images\IR\',num2str(z),'.png'];  
%=================================================================================================
  %% 读图
%=====================================================
img1 = double((imread(p1)));  
img2 = double((imread(p2)));

img1 = 255/(max(img1(:))-min(img1(:)))*(img1-min(img1(:))); 
img2 =255/(max(img2(:))-min(img2(:)))*(img2-min(img2(:)));
img1 = img1/255; img2 = img2/255;
figure;imshow(img1);
figure;imshow(img2);
%===========================================================
%% 定义变量及参数
%===============================================
nLevel = 4;
[m,n]=size(img1);
Sc1 = cell(1,nLevel+1);    Sc2 = cell(1,nLevel+1);%多尺度分解
MG1 = cell(1,nLevel+1);    MG2 = cell(1,nLevel+1);%增益函数
Lc1 = cell(1,nLevel);      Lc2 = cell(1,nLevel);%带通图像
L1 =  cell(1,nLevel);      L2 =  cell(1,nLevel);%局部带通对比度图像
TTA1 = cell(1,nLevel);     TTA2 = cell(1,nLevel);    TTA = cell(1,nLevel); % 感知对比度图像
DTA1 = cell(1,nLevel);     DTA2 = cell(1,nLevel);    DTA = cell(1,nLevel);% 去噪以及抑制过饱和后对比度图像
a = [1.40,1.15,1.04,1.15,1.35,1.93];
%=====================================

%% 多尺度变换
%=================================================================================
Sc1{1} = img1;  Sc2{1} = img2;
lambda1 =0.1;
Sc1{2} = SAEP(Sc1{1}, lambda1,4);Sc2{2} = SAEP(Sc2{1}, lambda1,4);
lambda2 =1;
Sc1{3} = SAEP(Sc1{2}, lambda2,8,4);Sc2{3} = SAEP(Sc2{2}, lambda2,8,4);
lambda3 =1.9;
Sc1{4} = SAEP(Sc1{3}, lambda3,16,8);Sc2{4} = SAEP(Sc2{3}, lambda3,16,8);
lambda4 =2.8;
Sc1{5} = SAEP(Sc1{4}, lambda4,32,16);Sc2{5} = SAEP(Sc2{4}, lambda4,32,16);

y1=0;y2=0;
for i = 1 : 4
    Lc1{i} = (255.*Sc1{i} - 255.*Sc1{i+1});Lc2{i} =(255.*Sc2{i} - 255.*Sc2{i+1});
    
%增益函数
    MG1{i} = gainv(255.*Sc1{i+1});  
    MG2{i} = gainr(255.*Sc2{i+1}); 

%得到适应的对比度信号 
    L1{i} = MG1{i}.*Lc1{i};  L2{i} =MG2{i}.*Lc2{i}; 

%% transducer处理及细节调整
%-------------------------------------
%非线性对比度转换
TTA1{i} = transducer(L1{i},a(i));   TTA2{i} = transducer(L2{i},a(i));

%视觉对比度调整，包括噪声和过曝光的抑制
DTA1{i} = detajf1(TTA1{i},Sc1{1},i);  
DTA2{i} = detajf1(TTA2{i},Sc2{1},i);

%% 融合系数准备
%--------------------------------
d1 = DTA1{i}; d2 = DTA2{i};

%*********Zhou**
% if i<2
%     d1 = guidedfilter(DTA1{i},DTA1{i},5,5); d2 = guidedfilter(DTA2{i},DTA2{i},5,5);
% end
%*********Zhou**
 
y1 = y1+d1; y2 = y2+d2;

end

% % ========================================================
%% 低通图像处理
%===================================================
LC1 = sbgain3(255.*Sc1{5},1);   LC2 = sbgain3(255.*Sc2{5},2);  %低通适应性处理
g1 = LC1+DTA1{4}; g2 = LC2+DTA2{4};
c=2.^(4+0);
ws1 = saliency_detail(g1,c*3); ws2 = saliency_detail(g2,c*3);
w1 = ws1; w2 = ws2;
w = w1./(w1+w2);
b = fspecial('gaussian',[6*c+1,6*c+1],c);
w = imfilter(w,b,'replicate');
LC = w.*LC1+(1-w).*LC2;

%% 带通参数融合
%==========================================================================
x1 = LC1+DTA1{nLevel};x2 = LC2+DTA2{nLevel};

for i = nLevel:-1:2
    if i<nLevel
        y1 = y1-DTA1{i+1}; y2 = y2-DTA2{i+1};  %像素级信息
        x1 = x1+DTA1{i};  x2 = x2+DTA2{i}; %结构显著性
    end
        
    DTA{i} = fusion_last(DTA1{i},DTA2{i},x1,x2,y1,y2,i,img2.*255,Sc2{i+1}.*255,Lc2{i});
end

DTA{1} = fusion_last(DTA1{1},DTA2{1},x1,x2,DTA1{1},DTA2{1},1,img2.*255,Sc2{2}.*255,Lc2{1});
% %========================================================================================================
%% 图像重构
%================================================
a = [1.40,1.15,1.04,1.15,1.35,1.93]; %非线性转换系数p
for i = nLevel:-1:1  
    %带通inv-transducer
    L{i} = in_transducer(DTA{i},a(i));

    %带通除以低通gain参数
    AL = L{i}./gainv(LC);
  
    %相加得到上一层次低通图像
    LC = AL + LC;

end

%================================


%% 图像显示
%==========================
img = LC; 
% img=ImRegular(LC);
img0 = uint8(img);
figure;imshow(img0);  

% imwrite(img/255,['.\results\',num2str(z),'.png']);

% % % ==========================================
 end
