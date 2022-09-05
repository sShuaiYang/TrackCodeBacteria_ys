function [type,Z,typeBacteria]=segementationImage(I,bacteriaNum)
sita=0:1:180;
bacteriaLength=65;
length=10;
width=5;
t=0;
[xSize,ySize]=size(I);
for i=1:size(sita,2)
    se=createSe(sita(i),length,width);
    I2=imerode(I,se);
    if max(max(I2))~=0
        cc=bwconncomp(I2);
        for j=1:cc.NumObjects
        t=t+1;
        erodePaper=false(xSize,ySize);
        erodePaper(cc.PixelIdxList{j})=1;
        erodePaper=imdilate(erodePaper,se);
        STATS=regionprops(erodePaper,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');
        ellipseInfo(t,:)=[STATS.Centroid(1)/bacteriaLength,STATS.Centroid(2)/bacteriaLength,STATS.MajorAxisLength/bacteriaLength,STATS.MinorAxisLength/bacteriaLength];
        orientationInfo(t,:)=STATS.Orientation/180*pi;
%         STATS=regionprops(erodePaper,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');
%         ellipseInfo(t,:)=[STATS.Centroid(1)/bacteriaLength,STATS.Centroid(2)/bacteriaLength];
        eachInfo(:,:,t)=erodePaper;
        end
    end
end
        D1=squareform(pdist(ellipseInfo,'euclidean'));
        sita=@(XI,XJ)((sin(XI-XJ)).^2);
        D2=squareform(pdist(orientationInfo, @(Xi,Xj)sita(Xi,Xj)));
Z=linkage(sqrt(D1.^2+D2),'ward');
c=cluster(Z,'maxclust',bacteriaNum);
for i=1:bacteriaNum
    finalPicture=false(xSize,ySize);
    numSequence=find(c==i);
    for j=1:numel(numSequence)
        type{i}(:,:,j)=eachInfo(:,:,numSequence(j));
        finalPicture=finalPicture|eachInfo(:,:,numSequence(j));
    end
    typeBacteria(:,:,i)=finalPicture;
end
imshow(imoverLap(im2uint8(I),typeBacteria))
end
function se=createSe(sita,length,width)
se=ones(width,length);
se=imrotate(se,sita);
end