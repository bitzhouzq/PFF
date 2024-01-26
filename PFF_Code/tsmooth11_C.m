
function S = tsmooth11_C(I,lambda,wtb)
    if (~exist('lambda','var'))
       lambda=0.01;
    end   
%     if (~exist('sigma','var'))
%        sigma=3.0;
%     end 
%     if (~exist('sharpness','var'))
%         sharpness = 1e-12;
%           sharpness = 0.02;
%    end
    if (~exist('maxIter','var'))
       maxIter=3;
    end    
%     I = im2double(I);
    x = I;
%     sigma_iter = sigma;
    lambda = lambda/2.0;
%     dec=1;
%     wtb = OriGrad_MD(I, sigma_iter, sigma0);
%     wtb = min(wtb, wtbi);
%     wtbo = wtb;
%     global bwx; global bwy;
%     bwx = ones(size(I)); bwy = ones(size(I));
    for iter = 1:maxIter
%         [wx, wy] = computeTextureWeights32(x, sigma_iter, sharpness);
%         [wx, wy, wxy] = computeTextureWeights33(x, sigma_iter, 8);
%         [wx, wy] = computeTextureWeights1(I, x, sigma_iter, sharpness);
%         [wx, wy] = computeTextureWeights34(x, sigma_iter, wtb);
        [wx, wy] = computeTextureWeights341(x, I, wtb, iter);
        
%         x = solveLinearEquation(I, wx, wy, lambda);
        x = FWLS_MM(I, x, wx, wy, lambda, 3, 4);

%         figure, imshow(x);
    end
    S = x;      
end

function [retx, rety] = computeTextureWeights341(fin, Iin, wtb, iter) 
   
% �ڶ��������
%    fin = log(fin+eps); 
%    Iin = log(Iin+eps); 
   
   fx = diff(fin,1,2);
   fx = padarray(fx, [0 1 0], 'post');
   fy = diff(fin,1,1);
   fy = padarray(fy, [1 0 0], 'post');
%     [fx,fy] = GradientMethod(fin, 'sobel');

   fxi = diff(Iin,1,2);
   fxi = padarray(fxi, [0 1 0], 'post');
   fyi = diff(Iin,1,1);
   fyi = padarray(fyi, [1 0 0], 'post');
%     [fxi,fyi] = GradientMethod(Iin, 'sobel');
   
%    iter = 10^(iter-1);
   global bwx; global bwy;
%    bwxt = ones(size(fx)); bwyt = ones(size(fy));
%    bfxi = abs(fxi) < abs(fx)-0.0000001; bfyi = abs(fyi) < abs(fy)-0.0000001;
%    bwxt(bfxi) = (abs(fx(bfxi)) - abs(fxi(bfxi)))./abs(fxi(bfxi));  bwyt(bfyi) = (abs(fy(bfyi)) - abs(fyi(bfyi)))./abs(fyi(bfyi)); 
%    bwxt(bfxi) = exp(-bwxt(bfxi)/(0.3));  bwyt(bfyi) = exp(-bwyt(bfyi)/(0.3));
% %    bwxt(bfxi) = 1;  bwyt(bfyi) = 1;
%    bwx = bwxt; bwy = bwyt;
   
   mgi = sqrt(fxi.^2+fyi.^2); mg = sqrt(fx.^2+fy.^2);
   bfxy = mgi < mg; %0.0005
   fx(bfxy) = fxi(bfxy); fy(bfxy) = fyi(bfxy);
   
   bwxt = ones(size(fx)); bwyt = ones(size(fy));
   bwxt(bfxy) = (mg(bfxy)-mgi(bfxy))./mgi(bfxy);
%    if iter>4
%        mu = 1000;
%    else
       mu = 0.1*(2^(iter-1)); %0.1*(iter-1) 
%    end
   bwxt(bfxy) = exp(-bwxt(bfxy)/mu); bwyt = bwxt;
%    bwxt(bfxy) = 1;  bwyt(bfxy) = 1;
   bwx = bwxt; bwy = bwyt;
%    figure, imshow(bwx)
  
%    dg = max(mg-mgi-0.0001,0);
%    b = mgi>0.0001; dg(b) = 0;
%    figure, imshow(dg*1e5);
      
   vareps = 2e-3; %1e-12 %2e-3
   
%    wtx = sum(abs(fx),3)/size(fin,3);
%    wty = sum(abs(fy),3)/size(fin,3);

%    h = fspecial('gaussian', [3, 3], 0.5);
%    fin = imfilter(fin, h, 'symmetric');

%    fin = log(fin+eps); 
%    fx = diff(fin,1,2);
%    fx = padarray(fx, [0 1 0], 'post');
%    fy = diff(fin,1,1);
%    fy = padarray(fy, [1 0 0], 'post');
   wtx = sum(sqrt(fx.^2+fy.^2),3)/size(fin,3);
   wty = wtx;
%    figure, imshow(wtx);
%    wtx = abs(fx); wty = abs(fy);

%    wtbx = max((wtx.*wtb).*bwx,vareps).^(-1);
%    wtby = max((wty.*wtb).*bwy,vareps).^(-1);   
   wtbx = max((wtx.*wtb),vareps).^(-1);
   wtby = max((wty.*wtb),vareps).^(-1);
   wtbx = wtbx.*((bwx+1e-6).^(-1));
   wtby = wtbx;
   
   retx = wtbx;
   rety = wtby;

   
   retx(:,end) = 0;
   rety(end,:) = 0;
end