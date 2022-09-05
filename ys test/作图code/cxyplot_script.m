figure,scatter(dataHist_sf.intensityCollect,dataHist_mS.intensityCollect,50,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.1)

temp=dataHist_mS.intensityCollect;
temp1=dataHist_sf.intensityCollect;

b=temp1<=1500 & temp>3 & temp<=300;%直接用逻辑矩阵 方便

% b=ones(numel(temp1),1);
% for i=1:numel(temp1)
%     if temp1(i)>2000|temp(i)<=0|temp(i)>300
%         b(i)=0;
%     end
% end
b=logical(b);

a=[];
a(1,:)=dataHist_sf.intensityCollect(b);
a(2,:)=dataHist_mS.intensityCollect(b);
a=a';
figure,scatter(a(:,1),a(:,2),50,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.1);
n=70;
[N C]=hist3(a,[n,n]);%二维平面格子里找到数据点的数目
N(find(N==0))=NaN;
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
    
% cmap=colormap(parula(213));
I=max(scdata(:,3));

cmap=ones(I,3);
cmap(end,:)=[0,0.45,0.74];%[0.74,0.74,0.74]
cmap(:,1)=0.74:-0.74/(I-1):0;
cmap(:,2)=0.74:-(0.74-0.45)/(I-1):0.45;
cmap(1:end-1,3)=cmap(1:end-1,3).*0.74;

clor=[];
for i=1:size(scdata,1)
    clor(i,:)=cmap(scdata(i,3),:);    
end
sz = 50;
figure,scatter(scdata(:,1),scdata(:,2),sz,clor,'filled')

figure,scatterhist(dataAll{1}(:,1),dataAll{1}(:,2),'Kernel','on')
figure,scatterhist(dataAll{1}(:,1),dataAll{1}(:,2),'Kernel','overlay')

for i=1:9
figure,scatterhist(dataAll{i}(:,1),dataAll{i}(:,2),'Kernel','on');ylabel('mScarlet');xlabel('sfGFP');
end

 rg=[];
for i=1:13
rg(:,i)=dataAll{i}(1:3000,2)./dataAll{i}(1:3000,1);
end

