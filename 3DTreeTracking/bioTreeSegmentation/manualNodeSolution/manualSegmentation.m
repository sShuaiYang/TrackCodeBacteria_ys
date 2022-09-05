function outResultImage=manualSegmentation(pixelIdxListIn,pixelIdxListOut,is2Node,imageSize)
disp([numel(pixelIdxListIn),numel(pixelIdxListOut)])
color=colormap(jet(20));
MIJ.run('Close All')
close all
xyMinMax=getAllXyMinMax(pixelIdxListIn,pixelIdxListOut,imageSize);
[inMask,~,h1]=getColorMarkImage(pixelIdxListIn,color,imageSize,xyMinMax);
movegui(h1,[1,-1]);
[outMask,xyMin,h2]=getColorMarkImage(pixelIdxListOut,color,imageSize,xyMinMax,is2Node);
movegui(h2,[-1,-1]);
[x,y]=find(inMask==1);
centroidIn=[mean(x),mean(y)];
[x,y]=find(outMask==1);
centroidOut=[mean(x),mean(y)];
dis=-centroidIn+centroidOut;
MIJ.createImage(im2uint8(outMask))
drawnow;commandwindow;
stateInput=input('which situation:');
switch stateInput
    case 'pass'
        clear outResultImage
        outResultImage.state='pass';
        return
    case 'sep'
        iNew=1;
        image{1}=[];
        outResultImage.state='sep';
        image{iNew}=input(strcat('Please Find The ',num2str(iNew),' Piece in ImageJ:___'));
        pause(0.1)
        MIJ.run('Close')
        outResultImage.pixelIdxList{iNew}=xy2Idx(xyMin,image{iNew},imageSize);
        while ~strcmp(image{iNew},'over')
            outResultImage.outNum(iNew)=input(strcat('Please Input The Matched outNum:___'));
            outResultImage.pixelIdxList{iNew}=xy2Idx(xyMin,image{iNew},imageSize);
            iNew=iNew+1;
            image{iNew}=input(strcat('Please Find The ',num2str(iNew),' Piece in ImageJ:___'));
            pause(0.1)
            MIJ.run('Close')
        end
        existNum=unique(outResultImage.outNum);
        for iOut=1:numel(pixelIdxListOut)
            if ~ismember(iOut,existNum)
                outResultImage.outNum(end+1)=iOut;
                outResultImage.pixelIdxList{end+1}=pixelIdxListOut{iOut};
            end
        end 
        return
    case 'mistake'
        clear outResultImage
        outResultImage.state='mistake';
        return
    case 'mustDiv'
        clear outResultImage
        outResultImage.state='mustDiv';
        return
    case 'deal'
        rightDeal=0;
        outResultImage.outNum=zeros(1,numel(pixelIdxListIn));
        inRemain=1:numel(pixelIdxListIn);
        outRemain=1:numel(pixelIdxListOut);
        inRemain1=inRemain;
        outRemain1=outRemain;
        iResult=1;
        outResultImage.matchResult{iResult}=input(strcat('Please Input The Matched pair [inNum,0,outNum] :___'));
        while rightDeal==0
            while ~strcmp(outResultImage.matchResult{iResult},'over')
                if strcmp(outResultImage.matchResult{iResult},'rest')
                    inSet=[];
                    outSet=[];
                    for i=1:numel(outResultImage.matchResult)-1
                        index=find(outResultImage.matchResult{i}==0);
                        inSet=[inSet,outResultImage.matchResult{i}(1:index-1)];
                        outSet=[outSet,outResultImage.matchResult{i}(index+1:end)];
                    end
                    diffIn=setdiff(1:numel(pixelIdxListIn),inSet);
                    diffOut=setdiff(1:numel(pixelIdxListOut),outSet);
                    outResultImage.matchResult{iResult}=[diffIn,0,diffOut];
                end
                if strcmp(outResultImage.matchResult{iResult},'reback')
                    iResult=iResult-1;
                    index=find(outResultImage.matchResult{iResult}==0);
                    inRemain1=[inRemain1,outResultImage.matchResult{iResult}(1:index-1)];
                    outRemain1=[outRemain1,outResultImage.matchResult{iResult}(index+1:end)];
                    MIJ.run('Close All')
                    close all
                    [~,h1]=getColorMarkImageNew(pixelIdxListIn,imageSize,xyMinMax,inRemain1,color);
                    movegui(h1,[1,-1]);
                    [~,h2]=getColorMarkImageNew(pixelIdxListOut,imageSize,xyMinMax,outRemain1,color);
                    movegui(h2,[-1,-1]);
                    drawnow;commandwindow;
                    outResultImage.matchResult{iResult}=input(strcat('Please Input The Matched pair [inNum,0,outNum] :___'));
                else
                    MIJ.run('Close All')
                    close all
                    index=find(outResultImage.matchResult{iResult}==0);
                    inRemain1=setdiff(inRemain1,outResultImage.matchResult{iResult}(1:index-1));
                    outRemain1=setdiff(outRemain1,outResultImage.matchResult{iResult}(index+1:end));
                    [~,h1]=getColorMarkImageNew(pixelIdxListIn,imageSize,xyMinMax,inRemain1,color);
                    movegui(h1,[1,-1]);
                    [~,h2]=getColorMarkImageNew(pixelIdxListOut,imageSize,xyMinMax,outRemain1,color);
                    movegui(h2,[-1,-1]);
                    drawnow;commandwindow;
                    iResult=iResult+1;
                    outResultImage.matchResult{iResult}=input(strcat('Please Input The Matched pair [inNum,0,outNum] :___'));
                end
            end
            outResultImage.matchResult(iResult)=[];
            rightDeal=1;
            inSet=[];
            outSet=[];
            for i=1:numel(outResultImage.matchResult)
                index=find(outResultImage.matchResult{i}==0);
                if i>=2 && (any(ismember(outResultImage.matchResult{i}(1:index-1),inSet)) || any(ismember(outResultImage.matchResult{i}(index+1:end),outSet)))
                    rightDeal=0;
                    break
                end
                inSet=[inSet,outResultImage.matchResult{i}(1:index-1)];
                outSet=[outSet,outResultImage.matchResult{i}(index+1:end)];
            end
            if any(~ismember(inSet,1:numel(pixelIdxListIn))) || any(~ismember(outSet,1:numel(pixelIdxListOut)))
                rightDeal=0;
            end
            if rightDeal==0
                disp('inputWrong    tryAgain')
                outResultImage=outResultImageOri;
                iResult=iResultOri;
                MIJ.run('Close All')
                close all
                [~,h1]=getColorMarkImageNew(pixelIdxListIn,imageSize,xyMinMax,inRemain,color);
                movegui(h1,[1,-1]);
                [~,h2]=getColorMarkImageNew(pixelIdxListOut,imageSize,xyMinMax,outRemain,color);
                movegui(h2,[-1,-1]);
                drawnow;commandwindow;
                inRemain1=inRemain;
                outRemain1=outRemain;
                outResultImage.matchResult{iResult}=input(strcat('Please Input The Matched pair [inNum,0,outNum] :___'));
            end
            outResultImage.state='deal';
        end
    case 'deal1'
        rightDeal=0;
        outResultImage.outNum=zeros(1,numel(pixelIdxListIn));
        [iResultOri,outResultImageOri,inRemain,outRemain]=autoMatchingProcessing(pixelIdxListIn,pixelIdxListOut,imageSize,xyMinMax,3,dis);
        inRemain1=inRemain;
        outRemain1=outRemain;
        iResult=iResultOri;
        outResultImage=outResultImageOri;
        outResultImage.matchResult{iResult}=input(strcat('Please Input The Matched pair [inNum,0,outNum] :___'));
        while rightDeal==0
            while ~strcmp(outResultImage.matchResult{iResult},'over')
                if strcmp(outResultImage.matchResult{iResult},'rest')
                    inSet=[];
                    outSet=[];
                    for i=1:numel(outResultImage.matchResult)-1
                        index=find(outResultImage.matchResult{i}==0);
                        inSet=[inSet,outResultImage.matchResult{i}(1:index-1)];
                        outSet=[outSet,outResultImage.matchResult{i}(index+1:end)];
                    end
                    diffIn=setdiff(1:numel(pixelIdxListIn),inSet);
                    diffOut=setdiff(1:numel(pixelIdxListOut),outSet);
                    outResultImage.matchResult{iResult}=[diffIn,0,diffOut];
                end
                if strcmp(outResultImage.matchResult{iResult},'reback')
                    iResult=iResult-1;
                    index=find(outResultImage.matchResult{iResult}==0);
                    inRemain1=[inRemain1,outResultImage.matchResult{iResult}(1:index-1)];
                    outRemain1=[outRemain1,outResultImage.matchResult{iResult}(index+1:end)];
                    MIJ.run('Close All')
                    close all
                    [~,h1]=getColorMarkImageNew(pixelIdxListIn,imageSize,xyMinMax,inRemain1,color);
                    movegui(h1,[1,-1]);
                    [~,h2]=getColorMarkImageNew(pixelIdxListOut,imageSize,xyMinMax,outRemain1,color);
                    movegui(h2,[-1,-1]);
                    drawnow;commandwindow;
                    outResultImage.matchResult{iResult}=input(strcat('Please Input The Matched pair [inNum,0,outNum] :___'));
                else
                    MIJ.run('Close All')
                    close all
                    index=find(outResultImage.matchResult{iResult}==0);
                    inRemain1=setdiff(inRemain1,outResultImage.matchResult{iResult}(1:index-1));
                    outRemain1=setdiff(outRemain1,outResultImage.matchResult{iResult}(index+1:end));
                    [~,h1]=getColorMarkImageNew(pixelIdxListIn,imageSize,xyMinMax,inRemain1,color);
                    movegui(h1,[1,-1]);
                    [~,h2]=getColorMarkImageNew(pixelIdxListOut,imageSize,xyMinMax,outRemain1,color);
                    movegui(h2,[-1,-1]);
                    drawnow;commandwindow;
                    iResult=iResult+1;
                    outResultImage.matchResult{iResult}=input(strcat('Please Input The Matched pair [inNum,0,outNum] :___'));
                end
            end
            outResultImage.matchResult(iResult)=[];
            rightDeal=1;
            inSet=[];
            outSet=[];
            for i=1:numel(outResultImage.matchResult)
                index=find(outResultImage.matchResult{i}==0);
                if i>=2 && (any(ismember(outResultImage.matchResult{i}(1:index-1),inSet)) || any(ismember(outResultImage.matchResult{i}(index+1:end),outSet)))
                    rightDeal=0;
                    break
                end
                inSet=[inSet,outResultImage.matchResult{i}(1:index-1)];
                outSet=[outSet,outResultImage.matchResult{i}(index+1:end)];               
            end
            if any(~ismember(inSet,1:numel(pixelIdxListIn))) || any(~ismember(outSet,1:numel(pixelIdxListOut)))
                rightDeal=0;
            end
            if rightDeal==0
                disp('inputWrong    tryAgain')
                outResultImage=outResultImageOri;
                iResult=iResultOri;
                MIJ.run('Close All')
                close all
                [~,h1]=getColorMarkImageNew(pixelIdxListIn,imageSize,xyMinMax,inRemain,color);
                movegui(h1,[1,-1]);
                [~,h2]=getColorMarkImageNew(pixelIdxListOut,imageSize,xyMinMax,outRemain,color);
                movegui(h2,[-1,-1]);
                drawnow;commandwindow;
                inRemain1=inRemain;
                outRemain1=outRemain;
                outResultImage.matchResult{iResult}=input(strcat('Please Input The Matched pair [inNum,0,outNum] :___'));
            end
            outResultImage.state='deal';
        end
    case 'deal2'
        rightDeal=0;
        outResultImage.outNum=zeros(1,numel(pixelIdxListIn));
        [iResultOri,outResultImageOri,inRemain,outRemain]=autoMatchingProcessing(pixelIdxListIn,pixelIdxListOut,imageSize,xyMinMax,8,dis);
        inRemain1=inRemain;
        outRemain1=outRemain;
        iResult=iResultOri;
        outResultImage=outResultImageOri;
        outResultImage.matchResult{iResult}=input(strcat('Please Input The Matched pair [inNum,0,outNum] :___'));
        while rightDeal==0
            while ~strcmp(outResultImage.matchResult{iResult},'over')
                if strcmp(outResultImage.matchResult{iResult},'rest')
                    inSet=[];
                    outSet=[];
                    for i=1:numel(outResultImage.matchResult)-1
                        index=find(outResultImage.matchResult{i}==0);
                        inSet=[inSet,outResultImage.matchResult{i}(1:index-1)];
                        outSet=[outSet,outResultImage.matchResult{i}(index+1:end)];
                    end
                    diffIn=setdiff(1:numel(pixelIdxListIn),inSet);
                    diffOut=setdiff(1:numel(pixelIdxListOut),outSet);
                    outResultImage.matchResult{iResult}=[diffIn,0,diffOut];
                end
                if strcmp(outResultImage.matchResult{iResult},'reback')
                    iResult=iResult-1;
                    index=find(outResultImage.matchResult{iResult}==0);
                    inRemain1=[inRemain1,outResultImage.matchResult{iResult}(1:index-1)];
                    outRemain1=[outRemain1,outResultImage.matchResult{iResult}(index+1:end)];
                    MIJ.run('Close All')
                    close all
                    [~,h1]=getColorMarkImageNew(pixelIdxListIn,imageSize,xyMinMax,inRemain1,color);
                    movegui(h1,[1,-1]);
                    [~,h2]=getColorMarkImageNew(pixelIdxListOut,imageSize,xyMinMax,outRemain1,color);
                    movegui(h2,[-1,-1]);
                    drawnow;commandwindow;
                    outResultImage.matchResult{iResult}=input(strcat('Please Input The Matched pair [inNum,0,outNum] :___'));
                else
                    MIJ.run('Close All')
                    close all
                    index=find(outResultImage.matchResult{iResult}==0);
                    inRemain1=setdiff(inRemain1,outResultImage.matchResult{iResult}(1:index-1));
                    outRemain1=setdiff(outRemain1,outResultImage.matchResult{iResult}(index+1:end));
                    [~,h1]=getColorMarkImageNew(pixelIdxListIn,imageSize,xyMinMax,inRemain1,color);
                    movegui(h1,[1,-1]);
                    [~,h2]=getColorMarkImageNew(pixelIdxListOut,imageSize,xyMinMax,outRemain1,color);
                    movegui(h2,[-1,-1]);
                    drawnow;commandwindow;
                    iResult=iResult+1;
                    outResultImage.matchResult{iResult}=input(strcat('Please Input The Matched pair [inNum,0,outNum] :___'));
                end
            end
            outResultImage.matchResult(iResult)=[];
            rightDeal=1;
            inSet=[];
            outSet=[];
            for i=1:numel(outResultImage.matchResult)
                index=find(outResultImage.matchResult{i}==0);
                if i>=2 && (any(ismember(outResultImage.matchResult{i}(1:index-1),inSet)) || any(ismember(outResultImage.matchResult{i}(index+1:end),outSet)))
                    rightDeal=0;
                    break
                end
                inSet=[inSet,outResultImage.matchResult{i}(1:index-1)];
                outSet=[outSet,outResultImage.matchResult{i}(index+1:end)];               
            end
            if any(~ismember(inSet,1:numel(pixelIdxListIn))) || any(~ismember(outSet,1:numel(pixelIdxListOut)))
                rightDeal=0;
            end
            if rightDeal==0
                disp('inputWrong    tryAgain')
                outResultImage=outResultImageOri;
                iResult=iResultOri;
                MIJ.run('Close All')
                close all
                [~,h1]=getColorMarkImageNew(pixelIdxListIn,imageSize,xyMinMax,inRemain,color);
                movegui(h1,[1,-1]);
                [~,h2]=getColorMarkImageNew(pixelIdxListOut,imageSize,xyMinMax,outRemain,color);
                movegui(h2,[-1,-1]);
                drawnow;commandwindow;
                inRemain1=inRemain;
                outRemain1=outRemain;
                outResultImage.matchResult{iResult}=input(strcat('Please Input The Matched pair [inNum,0,outNum] :___'));
            end
            outResultImage.state='deal';
        end
    case 1
        outResultImage.state=[];
end
end
function xyMinMax=getAllXyMinMax(in,out,imageSize)
pixelAll=[];
for i=1:numel(in)
    pixelAll=[pixelAll;in{i}];
end
for i=1:numel(out)
    pixelAll=[pixelAll;out{i}];
end
[~,~,xyMinMax]=idx2Xy(pixelAll,imageSize);
end
function [maskImage,xyMin,h]=getColorMarkImage(pixelIdxList,color,imageSize,xyMinMax,is2Node)
inNum=numel(pixelIdxList);
colorNum=fix(linspace(1,20,inNum));
image=uint8(false(xyMinMax(3)-xyMinMax(1)+3,xyMinMax(4)-xyMinMax(2)+3));
imageAll=cat(3,image,image,image);
maskImage=false(xyMinMax(3)-xyMinMax(1)+3,xyMinMax(4)-xyMinMax(2)+3);
for i=1:inNum
    [in{i}.xyMin,in{i}.BWImage]=idx2Xy(pixelIdxList{i},imageSize,xyMinMax);
    in{i}.BWImage=bwmorph(in{i}.BWImage,'remove');
    in{i}.Centroid=regionprops(in{i}.BWImage,'centroid');
    image1=imageAll(:,:,1);
    image2=imageAll(:,:,2);
    image3=imageAll(:,:,3);
    image1(in{i}.BWImage==1)=255*color(colorNum(i),1);
    image2(in{i}.BWImage==1)=255*color(colorNum(i),2);
    image3(in{i}.BWImage==1)=255*color(colorNum(i),3);
    imageAll=cat(3,image1,image2,image3);
    maskImage=maskImage | in{i}.BWImage;
end
figure;h=imshow(imageAll);
set(gcf,'Position',[240 195 694 627]);
for i=1:inNum
    if nargin==5
        if is2Node(i)==0
            text(in{i}.Centroid(1).Centroid(1),in{i}.Centroid(1).Centroid(2),num2str(i),'Color',[0,1,0],'FontSize',10);
            continue
        end
    end
%     text(in{i}.Centroid.Centroid(1),in{i}.Centroid.Centroid(2),num2str(i),'Color',color(colorNum(i),:),'FontSize',10);
    text(in{i}.Centroid(1).Centroid(1),in{i}.Centroid(1).Centroid(2),num2str(i),'Color',[1,1,1],'FontSize',10);
end
xyMin=in{1}.xyMin;
end
function [xyMin,BWImageGain,xyMinMax]=idx2Xy(pixelIdxList,pictureSize,xyMinMax)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=round(pixelIdxList-(yresult-1)*xSize);
if nargin==2
    xMin=min(xresult);
    xMax=max(xresult);
    yMin=min(yresult);
    yMax=max(yresult);
    xyMinMax=[xMin,yMin,xMax,yMax];
else
    xMin=xyMinMax(1);
    yMin=xyMinMax(2);
    xMax=xyMinMax(3);
    yMax=xyMinMax(4);
end
xyMin=[xMin,yMin];
yresult2=yresult-yMin+1;
xresult2=xresult-xMin+1;
BWImage=false(xMax-xMin+1,yMax-yMin+1);
pixelIdxList2=xresult2+(yresult2-1)*(xMax-xMin+1);
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
BWImageGain(2:end-1,2:end-1)=BWImage;
BWImageGain=imfill(BWImageGain,'holes');
end
function pixelIdxList=xy2Idx(xyMin,BWImageGain,pictureSize)
%this function is used to convert a small BWImage to its original pixelIdx
BWImage=BWImageGain(2:end-1,2:end-1);
% BWImage=bwmorph(BWImage,'remove');
pixelIdxListOri=find(BWImage==1);
if size(pixelIdxListOri,2)~=1
    pixelIdxListOri=pixelIdxListOri';
end
smallxSize=size(BWImage,1);
yresult=ceil(pixelIdxListOri/smallxSize);
xresult=pixelIdxListOri-(yresult-1)*smallxSize;
xresult2=xresult+xyMin(1)-1;
yresult2=yresult+xyMin(2)-1;
xSize=pictureSize(1);
pixelIdxList=xresult2+(yresult2-1)*xSize;
end
function [iResult,outResultImage,inRemain,outRemain]=autoMatchingProcessing(pixelIdxListIn,pixelIdxListOut,imageSize,xyMinMax,threShold,dis)
outResultImage=[];
inputNum=numel(pixelIdxListIn);
for i=1:inputNum
[regionNum,xyMin,regionImage]=findRegionNum(pixelIdxListIn{i},imageSize,xyMinMax);
if regionNum>1
    iResult=0;
    inRemain=[];
    outRemain=[];
    outResultImage=[];
    return
else
    in{i}.xyMin=xyMin;
    in{i}.BWImage=regionImage{1};
end
end
outputNum=size(pixelIdxListOut,2);
for i=1:outputNum
    [regionNum,xyMin,regionImage]=findRegionNum(pixelIdxListOut{i},imageSize,xyMinMax);
    if regionNum>1
        iResult=0;
        inRemain=[];
        outRemain=[];
        outResultImage=[];
        return
    else
        out{i}.xyMin=xyMin;
        out{i}.BWImage=regionImage{1};
end
end
if nargin==6
    se=translate(strel(1),round(dis));
end
inNum=numel(in);
for i=1:inNum
    if ~isempty(in{i}.BWImage)
        if nargin==6
            Info=regionprops(imdilate(in{i}.BWImage,se),'Centroid','MajorAxisLength','MinorAxisLength','Orientation');
        else
            Info=regionprops(in{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');
        end
        xyMin=in{i}.xyMin;
        inInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
    else
        inInfo(i,:)=zeros(1,4);
    end
end
outNum=numel(out);
for i=1:numel(out)
    if nargin==6
        Info=regionprops(out{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');
    else
        Info=regionprops(out{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');
    end
    xyMin=out{i}.xyMin;
    outInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
end
distanceMat=zeros(inNum,outNum);
distanceMat(1:inNum,1:outNum)=pdist2(inInfo(:,1:4),outInfo(:,1:4));
nMatch=0;
inSeries=[];
outSeries=[];

% 最近匹配原则      % 2014-10-3修正为最小匹配
% for iLine=1:inNum
%     if min(distanceMat(iLine,:))<=threShold
%         outOrder=find(distanceMat(iLine,:)==min(distanceMat(iLine,:)));
%         if numel(outOrder)>=2
%             continue
%         end
%         nMatch=nMatch+1;
%         inSeries=[inSeries;iLine];
%         outSeries=[outSeries;outOrder];
%         outResultImage.matchResult{nMatch}=[iLine,0,outOrder];
%     end
% end

% 最小匹配 % 2014-10-3
for iLine=1:inNum
    [a,b]=find(distanceMat==min(min(distanceMat)));
    if min(min(distanceMat))<=3;
        if numel(a)>=2
            continue
        end
        nMatch=nMatch+1;
        inSeries=[inSeries;a];
        outSeries=[outSeries;b];
        outResultImage.matchResult{nMatch}=[a,0,b];
        distanceMat(a,:)=100;
        distanceMat(:,b)=100;
    end
end

% 算出最近的匹配之后，由这些匹配的质心的偏移量作为台子缓慢移动的偏差进行修正
newNum=-1;
while numel(inSeries)~=newNum
    newNum=numel(inSeries);
    if newNum==0
        inToOut=[0,0];
    else
        inCentroidUse=inInfo(inSeries,1:2);
        outCentroidUse=outInfo(outSeries,1:2);
        inCen=[mean(inCentroidUse(:,1)),mean(inCentroidUse(:,2))];
        outCen=[mean(outCentroidUse(:,1)),mean(outCentroidUse(:,2))];
        inToOut=outCen-inCen;   %质心偏差
    end
end
inRemain=setdiff(1:inNum,inSeries);
outRemain=setdiff(1:outNum,outSeries);
inInfoNew=inInfo(inRemain,:);
outInfoNew=outInfo(outRemain,:);
outInfoNew(:,1)=outInfoNew(:,1)-inToOut(1);
outInfoNew(:,2)=outInfoNew(:,2)-inToOut(2);
distanceMat=zeros(size(inInfoNew,1),size(outInfoNew,1));
distanceMat=pdist2(inInfoNew(:,1:4),outInfoNew(:,1:4));
% 这一轮找到那些孤独的配对者，也就是有一个最小，且跟其他值差距很大（>4）
if numel(outRemain)>=2
    for iLine=1:numel(inRemain)
        iDis=distanceMat(iLine,:);
        iDisSort=sort(iDis);
        if iDisSort(2)-iDisSort(1)>4 && iDisSort(1)<=threShold;
            iDisOut=find(iDis==min(iDis));
            nMatch=nMatch+1;
            inSeries=[inSeries;inRemain(iLine)];
            outSeries=[outSeries;outRemain(iDisOut)];
            outResultImage.matchResult{nMatch}=[inRemain(iLine),0,outRemain(iDisOut)];
            distanceMat(iLine,:)=100;
            distanceMat(:,iDisOut)=100;
        end
    end
    for iLine=1:numel(inRemain)
        if min(min(distanceMat))==100
            continue
        end
        [a,b]=find(distanceMat==min(min(distanceMat)));
        if min(min(distanceMat))<=threShold;
            if numel(a)>=2
                continue
            end
            nMatch=nMatch+1;
            inSeries=[inSeries;inRemain(a)];
            outSeries=[outSeries;outRemain(b)];
            outResultImage.matchResult{nMatch}=[inRemain(a),0,outRemain(b)];
            distanceMat(a,:)=100;
            distanceMat(:,b)=100;
        end
    end
end
inRemain=setdiff(1:inNum,inSeries);
outRemain=setdiff(1:outNum,outSeries);
iResult=nMatch+1;
MIJ.run('Close All')
color=colormap(jet(20));
[~,h1]=getColorMarkImageNew(pixelIdxListIn,imageSize,xyMinMax,inRemain,color);
movegui(h1,[1,-1]);
[~,h2]=getColorMarkImageNew(pixelIdxListOut,imageSize,xyMinMax,outRemain,color);
movegui(h2,[-1,-1]);
drawnow;commandwindow;
end
function [regionNum,xyMin,regionImage]=findRegionNum(pixelIdxList,pictureSize,xyMinMax)
% find how many regions there are in a input/output
if isempty(pixelIdxList)
    regionNum=0;
    xyMin=zeros(1,2);
    regionImage{1}=[];
else
    if nargin==2
        [xyMin,BWImage]=idx2Xy(pixelIdxList,pictureSize);
    else
        [xyMin,BWImage]=idx2Xy(pixelIdxList,pictureSize,xyMinMax);
    end
    CC=bwconncomp(BWImage);
    regionNum=CC.NumObjects;
    for i=1:regionNum
        pixelIdxList2=CC.PixelIdxList{i};
        BW=false(CC.ImageSize);
        BW(pixelIdxList2)=1;
        regionImage{i}=BW;
    end
end
end
function [maskImage,h]=getColorMarkImageNew(pixelIdxList,imageSize,xyMinMax,order,color)
% color=colormap(jet(20));
inNum=numel(pixelIdxList);
colorNum=fix(linspace(1,20,inNum));
image=uint8(false(xyMinMax(3)-xyMinMax(1)+3,xyMinMax(4)-xyMinMax(2)+3));
imageAll=cat(3,image,image,image);
maskImage=false(xyMinMax(3)-xyMinMax(1)+3,xyMinMax(4)-xyMinMax(2)+3);
for i=order
    [in{i}.xyMin,in{i}.BWImage]=idx2Xy(pixelIdxList{i},imageSize,xyMinMax);
    in{i}.BWImage=bwmorph(in{i}.BWImage,'remove');
    in{i}.Centroid=regionprops(in{i}.BWImage,'centroid');
    image1=imageAll(:,:,1);
    image2=imageAll(:,:,2);
    image3=imageAll(:,:,3);
    image1(in{i}.BWImage==1)=255*color(colorNum(i),1);
    image2(in{i}.BWImage==1)=255*color(colorNum(i),2);
    image3(in{i}.BWImage==1)=255*color(colorNum(i),3);
    imageAll=cat(3,image1,image2,image3);
    maskImage=maskImage | in{i}.BWImage;
end
figure;h=imshow(imageAll);
set(gcf,'Position',[240 195 694 627]);
for i=order
    text(in{i}.Centroid(1).Centroid(1),in{i}.Centroid(1).Centroid(2),num2str(i),'Color',[1,1,1],'FontSize',14);
end
end
%% Used in & out bwImage as the input var
% function outResultImage=manualSegmentation(in,out,outAll)
% inNum=numel(in);
% outNum=numel(out);
% MIJ.run('Close All')
% for i=1:inNum
%     in{i}.Centroid=regionprops(in{i}.BWImage,'centroid');
%     in{i}.Centroid=in{i}.Centroid(end:-1:1);
%     in{i}.Centroid=in{i}.Centroid+in{i}.xyMin-[1,1];
%     MIJ.createImage(im2uint8(in{i}.BWImage))
%     MIJ.run('Rename...',strcat('title=In',num2str(i),'.tif'));
% end
% for i=1:outNum
%     MIJ.createImage(im2uint8(out{i}.BWImage))
%     MIJ.run('Rename...',strcat('title=Out',num2str(i),'.tif'));
% end
% outAll=getOutAll(out);
% MIJ.createImage(im2uint8(outAll.BWImage))
% for i=1:outNum
%     outResultImage{i}=input(strcat('Please Input The ',num2str(i),' Out:___'));
%     pause(0.1)
%     MIJ.run('Close')
% end
% end
% function outAll=getOutAll(out)
% pixelIdxList=[];
% for i=1:numel(out)
%     pixelIdxList=[pixelIdxList;xy2Idx(out{i}.xyMin,out{i}.BWImage,[2160,2560])];
% end
% [outAll.xyMin,outAll.BWImage]=idx2Xy(pixelIdxList,[2160,2560]);
% end
% function [xyMin,BWImageGain]=idx2Xy(pixelIdxList,pictureSize)
% % this function is used to convert linearIdx to a small BWImage
% xSize=pictureSize(1);
% yresult=ceil(pixelIdxList/xSize);
% xresult=round(pixelIdxList-(yresult-1)*xSize);
% xMin=min(xresult);
% xMax=max(xresult);
% yMin=min(yresult);
% yMax=max(yresult);
% yresult2=yresult-yMin+1;
% xresult2=xresult-xMin+1;
% BWImage=false(xMax-xMin+1,yMax-yMin+1);
% pixelIdxList2=xresult2+(yresult2-1)*(xMax-xMin+1);
% BWImage(pixelIdxList2)=true;
% BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
% xyMin=[xMin,yMin];
% BWImageGain(2:end-1,2:end-1)=BWImage;
% % BWImageGain=imfill(BWImageGain,'holes');
% end
% function pixelIdxList=xy2Idx(xyMin,BWImageGain,pictureSize)
% %this function is used to convert a small BWImage to its original pixelIdx
% BWImage=BWImageGain(2:end-1,2:end-1);
% % BWImage=bwmorph(BWImage,'remove');
% pixelIdxListOri=find(BWImage==1);
% if size(pixelIdxListOri,2)~=1
%     pixelIdxListOri=pixelIdxListOri';
% end
% smallxSize=size(BWImage,1);
% yresult=ceil(pixelIdxListOri/smallxSize);
% xresult=pixelIdxListOri-(yresult-1)*smallxSize;
% xresult2=xresult+xyMin(1)-1;
% yresult2=yresult+xyMin(2)-1;
% xSize=pictureSize(1);
% pixelIdxList=xresult2+(yresult2-1)*xSize;
% end
