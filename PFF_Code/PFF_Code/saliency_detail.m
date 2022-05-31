function y=saliency_detail(x,r)

%Óëbase²ãÒ»Ñù
[hei, wid] = size(x);
N = boxfilter(ones(hei, wid), r);
mean_xx = boxfilter(x.^2,r)./N;
mean_x = boxfilter(x,r)./N;
y = abs(mean_xx-2*x.*mean_x+x.^2);
b = fspecial('gaussian',[6*5+1,6*5+1],5);
y = imfilter(y,b,'replicate');

end