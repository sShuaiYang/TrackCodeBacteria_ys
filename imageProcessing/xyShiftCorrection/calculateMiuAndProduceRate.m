function [gamaAll,miuAll,fpProduceAll,timeAll]=calculateMiuAndProduceRate(bioTree,gfpImage,rfpImage)
% 板底扫描计算细菌的生长速率以及荧光蛋白的生成速率，由原始程序GRratioAnalysis稍加改动而成
iPiece=0;
% alphaG=0.4722*2.5;   % for EGFP
% alphaG=0.4722*5.7;   % for SFGFP
% alphaR=0.0534*20;
alphaG=3.6602;
alphaR=8.4102;
% ratio=0.40;
% E=1-alphaR/alphaG/ratio;
gMature=1;
rMature=log(2)/100;
pieceInfo=[];
for iframe=1:size(bioTree,2)
    for iRoot=1:size(bioTree{iframe}.root,2)
        if numel(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList)>=8
            iPiece=iPiece+1;
            pieceInfo{iPiece}.tracePixel=bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList;
            for iCell=1:numel(pieceInfo{iPiece}.tracePixel)
                pieceInfo{iPiece}.length(iCell,1)=bioTree{iframe}.root{iRoot}.traceInfo.measurment{iCell}.MajorAxisLength;
            end
            pieceInfo{iPiece}.traceRange=iframe:(iframe+numel(pieceInfo{iPiece}.tracePixel)-1);
        end
    end
    for iNode=1:size(bioTree{iframe}.node,2)
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            if numel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList)>=8
                iPiece=iPiece+1;
                pieceInfo{iPiece}.tracePixel=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList;
                for iCell=1:numel(pieceInfo{iPiece}.tracePixel)
                    pieceInfo{iPiece}.length(iCell,1)=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iCell}.MajorAxisLength;
                end
                pieceInfo{iPiece}.traceRange=iframe:(iframe+numel(pieceInfo{iPiece}.tracePixel)-1);
            end
        end
    end
end
if isempty(pieceInfo)   
    gamaAll=[];
    miuAll=[];
    fpProduceAll=[];
    return
else
%     gamaAll=zeros(1,numel(pieceInfo));
%     miuAll=zeros(1,numel(pieceInfo));
gamaAll=[];
miuAll=[];
timeAll=[];
    fpProduceAll=zeros(1,numel(pieceInfo));
end
for i=1:numel(pieceInfo)
    stackInfo=pieceInfo{i}.traceRange;
    if stackInfo(1)>=600
        p=1;
    end
    pixelIdxList=pieceInfo{i}.tracePixel;
    meanGFP=[];
    meanRFP=[];
    for iTrace=1:numel(pixelIdxList);
        pixelInfo=pixelIdxList{iTrace};
        gImage=gfpImage(:,:,stackInfo(iTrace));
        rImage=rfpImage(:,:,stackInfo(iTrace));
        meanGFP(iTrace,1)=mean(double(gImage(pixelInfo)));
        meanRFP(iTrace,1)=mean(double(rImage(pixelInfo)));
    end
    gama=getGrowthRate(pieceInfo{i});
    meanGFP=meanGFP./alphaG;
    meanRFP=meanRFP./alphaR;
%     meanGFP=meanGFP+meanRFP.*E;
%     meanGFP=smooth(meanGFP);
%     meanRFP=smooth(meanRFP);
%     hold on;plot(meanGFP,'g');plot(meanRFP,'r')
%     plot(meanRFP./meanGFP)
%     close all
    deltaFP=meanGFP-meanRFP;
    timeSeries=didaClock(pieceInfo{i}.traceRange);
    deltaT=diff(timeSeries);
    deltaT(end+1)=deltaT(end);
    deltaGFP=diff(meanGFP);
    deltaGFP(end+1)=deltaGFP(end);
    deltaRFP=diff(meanRFP);
    deltaRFP(end+1)=deltaRFP(end);
    
    % gama就是由菌长计算出的生长速率，另由红绿比计算一下作为参照，假设mG很大
    miu=(-deltaRFP./deltaT+rMature*deltaFP)./meanRFP;
%     gama修正
%     miu2=-((meanGFP-meanRFP)+deltaGFP./deltaT/gMature-deltaRFP./deltaT/rMature)./(meanGFP/gMature-meanRFP/rMature);
    
    % P为荧光蛋白的生产速率
    d_deltaFP=diff(deltaFP);
    d_deltaFP(end+1)=d_deltaFP(end);
    fpProduce=(miu+rMature).*deltaFP+d_deltaFP./deltaT;
    
%     plot(gama);hold on;plot(miu,'g');plot(miu2,'r');
%     plot(fpProduce)
%     close all
%     gamaAll(i)=mean(gama);
gamaAll=[gamaAll;mean(gama)];
miuAll=[miuAll;mean(miu)];
%     miuAll(i)=mean(miu);
    fpProduceAll(i)=mean(fpProduce);
    timeAll=[timeAll;mean(timeSeries)];
end
end
function timeSeriesNew=didaClock(timeSeries,treeSize)
treeSize=655;
timer=1:treeSize;
timeSequence=(timer*2)';
timeSeriesNew=timeSequence(timeSeries);
end
function gama=getGrowthRate(pieceInfo)
% smooth 生长曲线之后，求斜率
timeSeries=didaClock(pieceInfo.traceRange);
lengthInfo=pieceInfo.length;
[lengthInfo,timeSeries]=lengthInfoCheck(lengthInfo,timeSeries);
lengthInfo=smooth(lengthInfo,7);
if numel(lengthInfo)<=2
    gama=[];
    return
end
gama=prepareCurveData(timeSeries,lengthInfo);
if ~isempty(gama)
%     gama=ones(numel(lengthInfo),1)*gama;
%     return
end
gama(1,1)=(lengthInfo(2)-lengthInfo(1))/(timeSeries(2)-timeSeries(1))/(lengthInfo(1));
gama(numel(lengthInfo),1)=(lengthInfo(end)-lengthInfo(end-1))/(timeSeries(end)-timeSeries(end-1))/(lengthInfo(end));
gama(2:(end-1),1)=(lengthInfo(3:end)-lengthInfo(1:(end-2)))./(timeSeries(3:end)-timeSeries(1:(end-2)))./(lengthInfo(2:(end-1)));
end
function [xyMin,BWImageGain]=idx2Xy(pixelIdxList,pictureSize)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=round(pixelIdxList-(yresult-1)*xSize);
xMin=min(xresult);
xMax=max(xresult);
yMin=min(yresult);
yMax=max(yresult);
yresult2=yresult-yMin+1;
xresult2=xresult-xMin+1;
BWImage=false(xMax-xMin+1,yMax-yMin+1);
pixelIdxList2=round(xresult2+(yresult2-1)*(xMax-xMin+1));
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
BWImageGain(2:end-1,2:end-1)=BWImage;
BWImageGain=imfill(BWImageGain,'holes');
xyMin=[xMin,yMin];
end
function pixelIdxList=xy2Idx(xyMin,BWImageGain,pictureSize)
%this function is used to convert a small BWImage to its original pixelIdx
BWImage=BWImageGain(2:end-1,2:end-1);
pixelIdxListOri=find(BWImage==1);
smallxSize=size(BWImage,1);
yresult=ceil(pixelIdxListOri/smallxSize);
xresult=pixelIdxListOri-(yresult-1)*smallxSize;
xresult2=xresult+xyMin(1)-1;
yresult2=yresult+xyMin(2)-1;
xSize=pictureSize(1);
pixelIdxList=xresult2+(yresult2-1)*xSize;
end
%% Fit: 'untitled fit 1'.
function gama= prepareCurveData(xData,yData)
xData=xData-xData(1);
% Set up fittype and options.
ft = fittype( 'a*exp(b*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.StartPoint = [yData(1) 0];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts); 
% plot(xData,yData);hold on;plot(xData,fitresult.a*exp(fitresult.b*xData))
% close all
gama=fitresult.b;
end
function  [gfpImage,rfpImage]=backGroundSubstract(gfpImage,rfpImage)
%首先要对红色荧光和绿色荧光图像去除背景
% normalizeR=100;
% normalizeG=100;
% for i=1:size(gfpImage,3)
%     rImage=rfpImage(:,:,i);
%     gImage=gfpImage(:,:,i);
%     rIndex=rImage(:);
%     gIndex=gImage(:);
%     rIndex=sort(rIndex);
%     gIndex=sort(gIndex);
%     rBack=mean(rIndex(1:end/5));
%     gBack=mean(gIndex(1:end/5));
%     rImage=rImage/(rBack/normalizeR)-normalizeR;
%     gImage=gImage/(gBack/normalizeG)-normalizeG;
%     rImage=uint16(double(rImage)./backGroundR);
%     gImage=uint16(double(gImage)./backGroundG);
%     rfpImage(:,:,i)=rImage;
%     gfpImage(:,:,i)=gImage;
% end
for i=1:size(gfpImage,3)
    rImage=rfpImage(:,:,i);
    gImage=gfpImage(:,:,i);
    rIndex=rImage(:);
    gIndex=gImage(:);
    rIndex=sort(rIndex);
    gIndex=sort(gIndex);
    rBack=mean(rIndex(1:end*2/3));
    gBack=mean(gIndex(1:end*2/3));
    rfpImage(:,:,i)=rfpImage(:,:,i)-rBack;
    gfpImage(:,:,i)=gfpImage(:,:,i)-gBack;
end
end
function [lengthInfo,timeSeries]=lengthInfoCheck(lengthInfo,timeSeries)
% flow cell里常常会遇到的问题，细菌可能会站立行走，导致长度减少很多，需要把这些数据去除掉
i=2;
for n=2:numel(lengthInfo)
    if i<=numel(lengthInfo)
       if  lengthInfo(i)-lengthInfo(i-1)<0
           lengthInfo(i)=[];
           timeSeries(i)=[];
           i=i-1;
       end
    end
    i=i+1;
end
% plot(timeSeries,lengthInfo)
end