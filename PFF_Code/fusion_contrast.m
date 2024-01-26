function y = fusion_contrast(tv,tr,xv,xr,yv,yr,a,ir,ir_b)
% tv,tr--���ںϵĿɼ��������Աȶ�ͼ��
% xv,xr--�ṹ��������Ϣ
% yv,yr--���ؼ���Ϣ
% a--�ֽ����
% ir--����ͼ��
%ir_b--�����߶ȷֽ�ͼ��

%% ����Ȩ��ͼ����ϵ��

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

%% �����������Ժͽṹ������
P_V = abs(yv);  P_R = abs(yr); % ������������

r=[2,6,12,24,48,96];

S_V = gradient_f(xv,r(a)); S_R = gradient_f(xr,r(a)); %�ṹ������

wv = 1*P_V+5*S_V;  wr = p.*(1*P_R+5*S_R); %�����Ժϲ�

d = wv-wr;

W = 1./(1+exp((0.1*2.^(4-a)).*(0-d))); %����Ȩ��

c = a; 
b1 = fspecial('gaussian',[3*c+1,3*c+1],c/2);

W = imfilter(W,b1,'replicate');

%% �ں�

y = W.*tv+(1-W).*tr;  %��Ȩ�ں�


end