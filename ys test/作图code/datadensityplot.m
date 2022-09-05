function [scdata] = datadensityplot(x,y,colorindex)
%�ھ��󷽸���ͳ�Ƹ�����������Ŀ���ܶ�,�������� ��һ��xֵ���ڶ���yֵ��
% X must be a matrix with two columns.
% x,y ȡ����������ע���ڶ��������µı仯����,������ͼ�任��ԭxy��x=10.^x;y=10.^y;
% x=log10(x);
% y=log10(y);

X(:,1)=x;
X(:,2)=y;

% colorindex=1;% 1������ɫ 2 ����Ȼ�ɫ

ms=25;%marker size

n=100; %Ĭ��ȡ100*100���Ӽ��㣬�����Զ���

[N,C]=hist3(X,[n,n]);%��άƽ��������ҵ����ݵ����Ŀ
N(N==0)=NaN;
k=1;
scdata=[];
for ix=1:n
    for iy=1:n
        if ~isnan(N(ix,iy))
            scdata(k,1)=C{1}(ix);
            scdata(k,2)=C{2}(iy);
            scdata(k,3)=N(ix,iy);
            k=k+1;
        end
    end
end

% [xi,yi]=meshgrid(C{1},C{2}); 
% zi = griddata(x,y,N,xi,yi);

scx=scdata(:,1);
scy=scdata(:,2);
c=scdata(:,3);


if colorindex==1
    %��ɫ
    color1=[0,0.45,0.74];
    color2=[0.95,0.98,1];
end
if colorindex==2
    %�Ȼ�ɫ
    color1=[0.85,0.33,0.1];
    color2=[1,0.92,0.89];
end

I=64;
%color�ֳ����Եȷֳ�64��
mymap(:,1)=linspace(color2(1),color1(1),I);
mymap(:,2)=linspace(color2(2),color1(2),I);
mymap(:,3)=linspace(color2(3),color1(3),I);
map = colormap(mymap);
ind = fix((c-min(c))/(max(c)-min(c))*(size(map,1)-1))+1;
h = [];
h=scatter(scx,scy,ms,map(ind,:),'filled');


end

