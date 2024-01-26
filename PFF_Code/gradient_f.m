function y = gradient_f(x,a) %����ʦ���򷽷� �ο����ģ�Multi-scale weighted gradient-based fusion for multi-focus images
%[p,q]=gradient(x);
[p,q]=GradientMethod(x, 'zhou'); %���ݶ�
pq=p+1i*q;
[c1, c2] = EigDecBlock(pq, a);
wt = sqrt((sqrt(c1)+sqrt(c2)).^2 + 0.1*(sqrt(c1)-sqrt(c2)).^2);
% wt = squeeze(wt);
% b = fspecial('gaussian',[6*a+1,6*a+1],a);
% y = imfilter(wt,b,'replicate');
y=wt;
end