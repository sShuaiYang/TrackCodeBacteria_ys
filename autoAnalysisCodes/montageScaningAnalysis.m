function montageScaningAnalysis()
% 针对Jzy & Xag 分析XyScan montage的板底细菌HJD1而设计的批处理程序
% 理论上程序流程分一下步骤
% 1.montage存储的时间序列的图片是分离的，首先进行整合，变成扫描视野个数的时间序列，并进行背景矫正（用罗丹明与钙黄绿素矫正）
% 2.进行图像二值转换，细菌识别
% 3.对图像xy漂移进行矫正，并将结果的红绿图像进行保存
% 4.由maskImage进行制作bioTree，并保存
% 5.从荧光图片中获得对应的荧光信息并计算对应的数值
% 数据结构为  目标文件夹--allResult--imageGFP--imageGFP001,imageGFP002
 
dirFile=uigetdir();
nameList=dir(dirFile);
for iField=1:numel(nameList)-2
    clc
    %%
%     iImage=image(:,:,iField);
    disp(iField)
    tic;
    fieldName=strcat(dirFile,'\',nameList(iField+2).name);
    cd(fieldName)
    dirMaskImage=strcat(fieldName,'\','maskImage');
    dirCyOFPchannel=strcat(fieldName,'\','CyOFP');
    dirGFPchannel=strcat(fieldName,'\','GFP');
    dirBioTreeResult=strcat(fieldName,'\','bioTreeResult');
    dirMiuResult=strcat(fieldName,'\','miu&Presult');
    mkdir(dirMaskImage)
    mkdir(dirBioTreeResult)
    mkdir(dirMiuResult)
    
%     % load进入不同的channel
%     frameInfo=load([dirCyOFPchannel,'\frameInfo']);
%     frameInfo=frameInfo.frameInfo;
%     timeInterval=etime(frameInfo(3,1:6),frameInfo(2,1:6))/60;
%     nameList1=dir(dirCyOFPchannel);
%     disp('loading')
%     % load CyOFP
%     n=0;
%     for iSmallImage=1:numel(nameList1)-3
%         if strcmp(nameList1(iSmallImage+3).name(1:5),'image')
%             n=n+1;
%             imageStack=import_tiff_stack([dirCyOFPchannel,'\',nameList1(iSmallImage+3).name]);
%             maskImage=fluoImageProcessing(imageStack);
%             maskImage=imerode(maskImage,ones(1));
%             maskImage=bwmorph(maskImage,'open');
%             maskImages(:,:,2*n-1)=maskImage;
%             maskImages(:,:,2*n)=maskImage;
%         end
%     end
%     % load GFP
%     nameList2=dir(dirGFPchannel);
%     n=0;
%     for iSmallImage=1:numel(nameList2)-3
%         if strcmp(nameList2(iSmallImage+3).name(1:5),'image')
%             n=n+1;
%             imageGFP(:,:,2*n-1)=import_tiff_stack([dirGFPchannel,'\',nameList2(iSmallImage+3).name]);
%             imageGFP(:,:,2*n)=import_tiff_stack([dirGFPchannel,'\',nameList2(iSmallImage+3).name]);
%         end
%     end
%     imageGFP=imageGFP(:,:,1:24);
%     maskImages=maskImages(:,:,1:24);
% %     % xy漂移的矫正   
% %     [maskImage,imageGFP,imageRFP]=fluoImageCorrection(maskImage,imageGFP,imageRFP);
% %     
% %     % 将结果储存一下
%      save(strcat(dirMaskImage,'\1'),'maskImages','imageGFP')
%      load(strcat(dirMaskImage,'\1'))
%     
%     % 从maskImage制作bioTree
%     [~,bioTree]=treeTrackingJob([],maskImages,1,size(maskImages,3),0);
%     bioTree{1}.imageSize=size(maskImages(:,:,1));
%     bioTree=leafRootRefineBioTree(bioTree,0);
%     bioTree=bioTreeSizeReduction(bioTree,0);
%     
%     bioTree=autoTidyBioTree(bioTree);
%     bioTree=agaroseBottomAutoSeg(bioTree);
%     bioTree=type2NodeReduction(bioTree,2);
%     bioTree=type3NodeReduction(bioTree,30);
%     bioTree=type4NodeReduction(bioTree);
%     
%     bioTree=bioTreeDigHoleAndMeasure(bioTree,1,bioTree{1}.imageSize(1),bioTree{1}.imageSize(2));   % 实心的
%     dirSmallTree=strcat(dirBioTreeResult,'\','1');
%     save(dirSmallTree,'bioTree')
%     
%     load(strcat(dirMaskImage,'\1'))
%     load(dirSmallTree);
%     [bacInfo]=mainInfomationGet(bioTree,imageGFP,timeInterval);
    dirBacInfo=strcat(fieldName,'\','bacInfo');
%     mkdir(dirBacInfo)
%     save(strcat(dirBacInfo,'\','bacInfo'),'bacInfo')
    
    load(strcat(dirBacInfo,'\','bacInfo'));
    [gama,miu,fpProduce,GRratio,gfp,rfp,isDead]=calculateMiuAndProduceRate(bacInfo);
    save(strcat(dirMiuResult,'\','3_只计算第一个点_new_1点'),'gama','miu','fpProduce','GRratio','gfp','rfp','isDead')
    toc;
end
end
%% 图像二值化程序
%% fluoImage processing
function [maskImage,imageStack]=fluoImageProcessing(imageStack)
% minIntensity=300;
% imageStack=imageStack-200;
for iStack=1:size(imageStack,3)
    image=imageStack(:,:,iStack);
    pixelInfo=image(:);
    pixelInfo=sort(pixelInfo);
    maxPixel=max(pixelInfo);
    image=image-100;
    edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
    gaussianFilter1=fspecial('gaussian',[10, 10],3);
    gaussianFilter2=fspecial('gaussian',[10, 10],2);
    image=uint8(double(image)/double(maxPixel)*255);
    image=imfilter(image,gaussianFilter1); % use Guassian blur filter process
    image=imfilter(image,edgeFilter); %use edgeFilter process
    image=imfilter(image,gaussianFilter2);
     thre=20;
    image=im2bw(image,thre/255);
    image=imclearborder(image);
    image=bwareaopen(image,100,4);
    image=logical(image);
    image=imfill(image,'holes');
    maskImage(:,:,iStack)=image;
%     maskImage(:,:,iStack)=bwmorph(image,'close');
end
end

%% 从maskImage矫正maskImage，gfpImage,rfpImage的程序
function [maskImage,gfpImage,rfpImage] = fluoImageCorrection(maskImage,gfpImage,rfpImage)
% 针对荧光图像的时间序列的xy矫正
% 输入变量为二值图像，核心矫正函数已经优化过
imageSize=size(maskImage);
if min(imageSize(1:2))>=1500;
    maskImageNew=maskImage(600:1100,600:1100,:);
else
    maskImageNew=maskImage;
end
maskImageNew1(:,:,2:size(maskImageNew,3)+1)=maskImageNew;
bestPosition=zeros(size(maskImage,3),2);
parfor i=2:size(maskImage,3)
    bestPosition(i,:)=caculateCrossCorrelationForImage(maskImageNew1(:,:,i),maskImageNew(:,:,i),15);
end
maskImage=imageCorrectionWithBestPosition(maskImage,bestPosition);
for i=1:size(maskImage,3)
    maskImage(:,:,i)=imclearborder(maskImage(:,:,i));
end
gfpImage=imageCorrectionWithBestPosition(gfpImage,bestPosition);
rfpImage=imageCorrectionWithBestPosition(rfpImage,bestPosition);
end
function bestPosition=caculateCrossCorrelationForImage(image1,image2,step)
% for calculater the cross correlation of two image
% step is the searching range
% backGround should be calculater or set by user
[x,y]=meshgrid((-step:step)',(-step:step)');
x=x(:);
y=y(:);
correlationMatrix=zeros(size(x,1),1);
parfor i=1:numel(x)
    se=translate(strel(1),[x(i),y(i)]);
    image2New=imdilate(image2,se);  % 巧用imdilate实现平移
    sumImage=image1 & image2New;    % 利用逻辑矩阵的乘法相当于&
    correlationMatrix(i)=sum(sum(sumImage)); % 直接对逻辑矩阵求和速度比较快
end
bestPosition=[x(correlationMatrix==max(correlationMatrix)),y(correlationMatrix==max(correlationMatrix))];
bestPosition=bestPosition(1,:);
end
function image=imageCorrectionWithBestPosition(image,bestPosition)
% 已知漂移量后进行的较正
bestPositionAccumulation=zeros(size(bestPosition));
for i=2:size(bestPosition,1)
    bestPositionAccumulation(i,:)=bestPosition(i,:)+bestPositionAccumulation(i-1,:);
end
for i=2:size(bestPosition,1)
    se=translate(strel(1),bestPositionAccumulation(i,:));
    image(:,:,i)=imdilate(image(:,:,i),se);
end
end

%% 制作bioTree过程中的子程序
function bioTree=bioTreeSizeReduction(bioTree,frameShift)
for iframe=frameShift+1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.is2Node==true
                bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList(end)=[];
%                 bioTree{iframe}.root{iroot}.traceInfo.measurment(end)=[];
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
                    bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList(end)=[];
%                     bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment(end)=[];
                end
            end
        end
    end
end
end

%% bioTreeMeasure 填充的
function bioTree=bioTreeDigHoleAndMeasure(bioTree,beginFrame,xSize,ySize)
parfor iframe=beginFrame:size(bioTree,2)
   if ~isempty(bioTree{iframe}.root)
       bioTree{iframe}.root=DigRoot(bioTree{iframe}.root,xSize,ySize);
   end
   if ~isempty(bioTree{iframe}.node)
       bioTree{iframe}.node=DigNode(bioTree{iframe}.node,xSize,ySize);
   end
   if ~isempty(bioTree{iframe}.leavies)
       bioTree{iframe}.leavies=DigLeaf(bioTree{iframe}.leavies,xSize,ySize);
   end
end
end
function bioRoot=DigRoot(bioRoot,xSize,ySize)
imageSize=[xSize,ySize];
for iRoot=1:size(bioRoot,2)
    pixelIdxList=bioRoot{iRoot}.rootPixelDetail;
    [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
    bioRoot{iRoot}.rootPixelDetail=pixelIdxList;
    bioRoot{iRoot}.rootMeasurment=pos;
    bioRoot{iRoot}.traceInfo.measurment=[];
    for iTrace=1:size(bioRoot{iRoot}.traceInfo.pixelIdxList,2)
        pixelIdxList=bioRoot{iRoot}.traceInfo.pixelIdxList{iTrace};
        [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
        bioRoot{iRoot}.traceInfo.measurment{iTrace}=pos;
        bioRoot{iRoot}.traceInfo.pixelIdxList{iTrace}=pixelIdxList;
    end
end
end
function bioLeaf=DigLeaf(bioLeaf,xSize,ySize)
imageSize=[xSize,ySize];
for iLeaf=1:size(bioLeaf,2)
    pixelIdxList=bioLeaf{iLeaf}.leaviesPixelDetail;
    [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
    bioLeaf{iLeaf}.leafPixelDetail=pixelIdxList;
    bioLeaf{iLeaf}.leafMeasurment=pos;
end
end
function bioNode=DigNode(bioNode,xSize,ySize)
imageSize=[xSize,ySize];
for iNode=1:size(bioNode,2)
    for iNodeOut=1:size(bioNode{iNode}.Out,2)
        bioNode{iNode}.Out{iNodeOut}.traceInfo.measurment=[];
        for iNodeTrace=1:size(bioNode{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList,2)
            pixelIdxList=bioNode{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace};
            [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
            bioNode{iNode}.Out{iNodeOut}.traceInfo.measurment{iNodeTrace}=pos;
            bioNode{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace}=pixelIdxList;
        end
    end
end
end
function [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize)
[xyMin,BWImage]=idx2XyMeasure(pixelIdxList,imageSize);
pos=regionprops(BWImage,'FilledArea','Centroid','Eccentricity','MajorAxisLength','Orientation','MinorAxisLength');
BWImage=imfill(BWImage,'holes');
for i=1:size(pos,1)
pos(i).Centroid(1)=pos(i).Centroid(1)+xyMin(2)-1;
pos(i).Centroid(2)=pos(i).Centroid(2)+xyMin(1)-1;
end
pixelIdxList=xy2IdxMeasure(xyMin,BWImage,imageSize);
end
function [xyMin,BWImageGain]=idx2XyMeasure(pixelIdxList,pictureSize)
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
% BWImageGain=bwmorph(BWImageGain,'close');
% BWImageGain=bwmorph(BWImageGain,'remove');
xyMin=[xMin,yMin];
end
function pixelIdxList=xy2IdxMeasure(xyMin,BWImageGain,pictureSize)
%this function is used to convert a small BWImage to its original pixelIdx
BWImage=BWImageGain(2:end-1,2:end-1);
% BWImage=bwmorph(BWImage,'remove');
pixelIdxListOri=find(BWImage==1);
smallxSize=size(BWImage,1);
yresult=ceil(pixelIdxListOri/smallxSize);
xresult=pixelIdxListOri-(yresult-1)*smallxSize;
xresult2=xresult+xyMin(1)-1;
yresult2=yresult+xyMin(2)-1;
xSize=pictureSize(1);
pixelIdxList=xresult2+(yresult2-1)*xSize;
end

%% bioTreeSegmentation
function bioTree=bioTreeSegmentation(bioTree)
bioTree=autoTidyBioTree(bioTree);
nodeSortPre=[];
[nodeSort,~]=findNode(bioTree);
while ~isequal(nodeSortPre,nodeSort)
    nodeSortPre=nodeSort;
    bioTree=type1NodeReduction(bioTree,1,25);
    bioTree=type2NodeReduction(bioTree,1);
    bioTree=mixMatchReducion(bioTree,5);
    bioTree=type3NodeReduction(bioTree,15);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
    [nodeSort,~]=findNode(bioTree);
end
while ~isequal(nodeSortPre,nodeSort)
    nodeSortPre=nodeSort;
    bioTree=type_InNodeReduction(bioTree,1);
    bioTree=type2NodeReduction(bioTree,1);
    bioTree=mixMatchReducion(bioTree,5);
    bioTree=type3NodeReduction(bioTree,15);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
    [nodeSort,~]=findNode(bioTree);
end
nodeSortPre=[];
while ~isequal(nodeSortPre,nodeSort)
    nodeSortPre=nodeSort;
    bioTree=type1NodeReduction(bioTree,2,25);
    bioTree=type_InNodeReduction(bioTree,2);
    bioTree=type2NodeReduction(bioTree,1);
    bioTree=mixMatchReducion(bioTree,5);
    bioTree=type3NodeReduction(bioTree,15);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
    [nodeSort,~]=findNode(bioTree);
end
nodeSortPre=[];
while ~isequal(nodeSortPre,nodeSort)
    nodeSortPre=nodeSort;
    bioTree=type1NodeReduction(bioTree,3,40);
    bioTree=type_InNodeReduction(bioTree,2);
    bioTree=type2NodeReduction(bioTree,1);
    bioTree=mixMatchReducion(bioTree,5);
    bioTree=type3NodeReduction(bioTree,15);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
    [nodeSort,~]=findNode(bioTree);
end
nodeSortPre=[];
while ~isequal(nodeSortPre,nodeSort)
    nodeSortPre=nodeSort;
    bioTree=type1NodeReduction(bioTree,3,40);
    bioTree=type_InNodeReduction(bioTree,2);
    bioTree=type2NodeReduction(bioTree,2);
    bioTree=mixMatchReducion(bioTree,10);
    bioTree=type3NodeReduction(bioTree,30);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
    [nodeSort,~]=findNode(bioTree);
end
bioTree=mergeNode(bioTree);
nodeSortPre=[];
while ~isequal(nodeSortPre,nodeSort)
    nodeSortPre=nodeSort;
    bioTree=type1NodeReduction(bioTree,3,40);
    bioTree=type_InNodeReduction(bioTree,2);
    bioTree=type2NodeReduction(bioTree,2);
    bioTree=mixMatchReducion(bioTree,10);
    bioTree=type3NodeReduction(bioTree,30);
    bioTree=type4NodeReduction(bioTree);
    [nodeSort,~]=findNode(bioTree);
end
bioTree=leafRootRefineBioTree(bioTree,0);
end


%% 提取计算时需要信息的主程序
function [bacInfo]=mainInfomationGet(bioTree,gfpImage,timeInterval)
% 提取计算时必要的信息，细菌的原始平均绿光，平均红光以及长度
% 板底扫描计算细菌的生长速率以及荧光蛋白的生成速率，由原始程序GRratioAnalysis稍加改动而成
iPiece=0;
pieceInfo=[];
for iframe=1:size(bioTree,2)
    for iRoot=1:size(bioTree{iframe}.root,2)
        if numel(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList)>=2
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
            if numel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList)>=2
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
bacInfo=[];
if isempty(pieceInfo)   
    return
end
for i=1:numel(pieceInfo)
    bacInfo{i}.lengthInfo=[];
    bacInfo{i}.meanGFP=[];
    bacInfo{i}.meanRFP=[];
    bacInfo{i}.time=[];
end
parfor i=1:numel(pieceInfo)
    stackInfo=pieceInfo{i}.traceRange;
    pixelIdxList=pieceInfo{i}.tracePixel;
    meanGFP=[];
    for iTrace=1:numel(pixelIdxList);
        pixelInfo=pixelIdxList{iTrace};
%         [xyMin,BWImageGain]=idx2Xy(pixelInfo,[512,512]);
%         pixelInfo=xy2Idx(xyMin,BWImageGain,[512,512]);
        gImage=gfpImage(:,:,stackInfo(iTrace));
        meanGFP(iTrace,1)=mean(double(gImage(pixelInfo))); 
    end
    [~,index]=find(mod(pieceInfo{i}.traceRange,2)==0);
    bacInfo{i}.lengthInfo=pieceInfo{i}.length(index);
    bacInfo{i}.meanGFP=meanGFP(index);
    bacInfo{i}.time=timeInterval*pieceInfo{i}.traceRange(index)/2;
    bacInfo{i}.time=bacInfo{i}.time';
    bacInfo{i}.pixelInfo=pixelIdxList{end}(index);
end
end

%% 计算各参数及最终数据的主程序
function [gamaAll,miuAll,fpProduceAll,GRratioAll,gfpAll,rfpAll,isDead]=calculateMiuAndProduceRate(bacInfo)
% 板底扫描计算细菌的生长速率以及荧光蛋白的生成速率，由原始程序GRratioAnalysis稍加改动而成
alphaG=0.3320*0.6;   % 2016-3-8红绿修正
alphaR=1.155*0.6;   % 2 result
gMature=log(2)/5;
rMature=log(2)/100;
if isempty(bacInfo)   
    gamaAll=[];
    miuAll=[];
    fpProduceAll=[];
    GRratioAll=[];
    gfpAll=[];
    rfpAll=[];
    isDead=[];
    timeAll=[];
    return
else
    gamaAll=[];
    miuAll=[];
    fpProduceAll=[];
    GRratioAll=[];
    gfpAll=[];
    rfpAll=[];
    isDead=[];
    timeAll=[];
end
for i=1:numel(bacInfo) 
% changed by jzy 2016.3.8 用第一二张的图片和生长速率来计算
    lengthInfo=bacInfo{i}.lengthInfo;
    timePoint=1;
    if numel(lengthInfo)<timePoint+1
        continue
    end
    meanGFP=bacInfo{i}.meanGFP*6510/6909;    %9941
    meanRFP=bacInfo{i}.meanRFP*11918/11227;   % 2016.03.20
    timeSeries=bacInfo{i}.time;
    gama=getGrowthRate(timeSeries,lengthInfo,timePoint);
    gama=gama';
    if isempty(gama)
        continue
    end
    meanGFP=meanGFP./alphaG/3;
    deltaT=diff(timeSeries);
    deltaGFP=meanGFP(2:end)-meanGFP(1:end-1); 
    % P为荧光蛋白的生产速率
    fpProduce=gama.*meanGFP(1:end-1)+deltaGFP./deltaT;
    isDead=0;
    gamaAll=[gamaAll;gama];
    fpProduceAll=[fpProduceAll;fpProduce];
    gfpAll=[gfpAll;meanGFP(1:end-1)];
    timeAll=[timeAll;timeSeries(1:end-1)];
end
end
function gama=getGrowthRate(timeSeries,lengthInfo,timePoint)
% smooth 生长曲线之后，求斜率
% lengthInfo=smooth(lengthInfo,7);
timeSeries=timeSeries;
lenT=numel(lengthInfo);
if numel(lengthInfo)<=1
    gama=[];
    return
end
% gama=prepareCurveData(timeSeries,lengthInfo);
% if gama<=0
%     gama=[];
%     return
% end
% gama=ones(numel(lenT),1)*gama;
% return
% gama(1,1)=(lengthInfo(2)-lengthInfo(1))/(timeSeries(2)-timeSeries(1))/(lengthInfo(1));
% while gama<0
%     i=i+1;
%     if i>=numel(lengthInfo)
%         break
%     end
for i=1:numel(lengthInfo)-1
    gama(i)=log2(lengthInfo(i+1)/lengthInfo(i))/(timeSeries(i+1)-timeSeries(i));
end
% gama(numel(lengthInfo),1)=(lengthInfo(end)-lengthInfo(end-1))/(timeSeries(end)-timeSeries(end-1))/(lengthInfo(end));
% gama(2:(end-1),1)=(lengthInfo(3:end)-lengthInfo(1:(end-2)))./(timeSeries(3:end)-timeSeries(1:(end-2)))./(lengthInfo(2:(end-1)));
end
%% Fit: 'untitled fit 1'.
function gama= prepareCurveData(xData,yData)
xData=xData-xData(1);
% Set up fittype and options.
yData1=log(yData);
ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.StartPoint = [0 0];
% Fit model to data.
% ft = fittype( 'a*exp(b*x)', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( ft );
% opts.StartPoint = [yData(1) 0];
[fitresult, gof] = fit( xData, yData1, ft, opts); 
% plot(xData,yData);hold on;plot(xData,fitresult.a*exp(fitresult.b*xData))
% close all

if gof.rsquare<=0.99
    gama=0;
    return
end
% fitresult
% plot(xData,yData1)
% close
gama=fitresult.a;
end
function xFP=sfGFPbleachCorrection(xFP)
r=0.002;
deltaFP=xFP*r;
for i=2:numel(deltaFP)-1
    deltaFP(i)=deltaFP(i)+deltaFP(i-1);
end
deltaFP=[0;deltaFP];
deltaFP(end)=[];
xFP=xFP+deltaFP;
end
function xFP=rFPbleachCorrection(xFP)
r=0.0139;
deltaFP=xFP*r;
for i=2:numel(deltaFP)-1
    deltaFP(i)=deltaFP(i)+deltaFP(i-1);
end
deltaFP=[0;deltaFP];
deltaFP(end)=[];
xFP=xFP+deltaFP;
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
pixelIdxList2=xresult2+(yresult2-1)*(xMax-xMin+1);
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
xyMin=[xMin,yMin];
BWImageGain(2:end-1,2:end-1)=BWImage;
BWImageGain=bwmorph(BWImageGain,'thin',inf);
BWImageGain=imdilate(BWImageGain,ones(3));
% BWImageGain=imfill(BWImageGain,'holes');
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
function manualTidyForSingleImage(dirFile)
fileList=dir(dirFile);
dirGFP=[dirFile,'\backGroundGFP.tif'];
dirRFP=[dirFile,'\backGroundRFP.tif'];
for i=1:numel(fileList)-2
    if ~strcmp(fileList(i+2).name(end-3:end),'.tif')
        cd([dirFile,'\',fileList(i+2).name])
        mkdir('original data')
        mkdir('information')
        movefile([dirFile,'\',fileList(i+2).name,'\*tif'],[dirFile,'\',fileList(i+2).name,'\original data']);
        nameList=dir();
        for iList=1:numel(nameList)-2
            if strcmp(nameList(iList+2).name(end-3:end),'.txt')
                movefile([dirFile,'\',fileList(i+2).name,'\',nameList(iList+2).name],[dirFile,'\',fileList(i+2).name,'\information\montageInfo.txt']);
            end
        end
        copyfile(dirGFP,[dirFile,'\',fileList(i+2).name,'\backGroundGFP.tif'])
        copyfile(dirRFP,[dirFile,'\',fileList(i+2).name,'\backGroundRFP.tif'])
    end
end
delete(dirGFP)
delete(dirRFP)
end
%% flowcell tracking auto segementation

%%% type1 node reduction
function [bioTree,usefulNode]=type1NodeReduction(bioTree,strength,finalOneToOneMatchThreShold)
% ##defalt finalOneToOneMatchThreShold=25 could be reach 40(in ASolution l120)
nodeCount=0;
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        nodeList=type1NodeinFrame(bioTree{iframe},iframe,bioTree);
        if ~isempty(nodeList)
            for iList=1:size(nodeList,1)
                nodeInfoIn=getInputInfo(bioTree,nodeList(iList,:));
                [traceInfo,usefulNode]=fullNodeType1Tracking(bioTree,nodeList(iList,:),strength,finalOneToOneMatchThreShold);
                nodeList(iList,4)=usefulNode;
                if usefulNode==0;
                    bioTree=type1NodeLinker(bioTree,nodeInfoIn,traceInfo);
                end
                nodeCount=nodeCount+1;
            end
            cantDivideNode=nodeList(:,4);
            nodeList(cantDivideNode==1,:)=[];
            if ~isempty(nodeList)
                bioTree=removeType1Node(bioTree,nodeList);
            end
        end
    end
end
end
function nodeList=type1NodeinFrame(bioTreeFrame,iframe,bioTree)
nodeList=[];
for iNode=1:size(bioTreeFrame.node,2)
    if size(bioTreeFrame.node{iNode}.In,2)>=1 && size(bioTreeFrame.node{iNode}.Out,2)==1
        %         if bioTreeFrame.node{iNode}.branchIndex~=12
        %             break
        %         end
        nodeInfo=[iframe,iNode];
        usefulNode=0;
        pixelIdxListIn=getInputMask(bioTree,nodeInfo);
        for iIn=1:size(pixelIdxListIn,2)
            [~,BWImage]=idx2Xy(pixelIdxListIn{iIn},bioTree{1}.imageSize);
            CC=bwconncomp(BWImage);
            regionNum=CC.NumObjects;
            if regionNum>=2
                usefulNode=1;
                break
            end
        end
        if usefulNode==0
            nodeList=[nodeList;[iframe,iNode]];
        end
    end
end
end
function bioTree=type1NodeLinker(bioTree,nodeInfoIn,traceInfo)
newNodeIn=0;
for iIn=1:size(nodeInfoIn,2)
    if isempty(traceInfo{iIn}.isBreakNode)
        if nodeInfoIn{iIn}.isNode==false
            rootInfo=nodeInfoIn{iIn}.rootInfo;
            if traceInfo{iIn}.is2Node==false;
                leafInfo=traceInfo{iIn}.leafInfo;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node=false;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[];
                bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                %                 leafIndex=size(bioTree{leafInfo(1)}.leavies,2);
                %                 leafInfo(2)=leafIndex+1;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=leafInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=false;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=rootInfo;
                if isempty(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList)
                    bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{1}=bioTree{rootInfo(1)}.root{rootInfo(2)}.rootPixelDetail;
                end
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{end};
            end
            if traceInfo{iIn}.is2Node==true;
                nodeInfoT=traceInfo{iIn}.nodeInfo;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=nodeInfoT;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=false;
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=[];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.rootInfo=rootInfo;
            end
        end
        if nodeInfoIn{iIn}.isNode==true
            nodeInfo=nodeInfoIn{iIn}.nodeInfo;
            if traceInfo{iIn}.is2Node==false;
                leafInfo=traceInfo{iIn}.leafInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node=false;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[];
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                %                 leafIndex=size(bioTree{leafInfo(1)}.leavies,2);
                %                 leafInfo(2)=leafIndex+1;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=leafInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=true;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=nodeInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList{end};
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo= leafInfo;
            end
            if traceInfo{iIn}.is2Node==true;
                nodeInfoT=traceInfo{iIn}.nodeInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo= nodeInfoT;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=true;
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=nodeInfo;
            end
        end
    end
    if traceInfo{iIn}.isBreakNode==1
        newNodeIn=newNodeIn+1;
        if newNodeIn==1
            breakNodeInfo=[traceInfo{iIn}.breakNodeInfo(1),size(bioTree{traceInfo{iIn}.breakNodeInfo(1)}.node,2)+1];
        end
        if nodeInfoIn{iIn}.isNode==false
            rootInfo=nodeInfoIn{iIn}.rootInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node=true;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[breakNodeInfo,newNodeIn];
            bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=[];
            bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.In{newNodeIn}.isNode=false;
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.In{newNodeIn}.rootInfo=rootInfo;
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.state=[];
        else
            if nodeInfoIn{iIn}.isNode==true
                preNodeInfo=nodeInfoIn{iIn}.nodeInfo;
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.is2Node=true;
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.nodeInfo=[breakNodeInfo,newNodeIn];
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.leafInfo=[];
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.In{newNodeIn}.isNode=1;
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.In{newNodeIn}.nodeInfo=preNodeInfo;
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.state=[];
            end
        end
        if  traceInfo{iIn}.is2Node==false
            leafInfo=traceInfo{iIn}.leafInfo;
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.is2Node=false;
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.leafInfo=leafInfo;
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.traceInfo.pixelIdxList=traceInfo{iIn}.afterBreakTrace.pixelIdxList;
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.nodeInfo=[];
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.state=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=true;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[breakNodeInfo,1];
        else
            if traceInfo{iIn}.is2Node==true
                nodeInfo=traceInfo{iIn}.nodeInfo;
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.is2Node=true;
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.nodeInfo=nodeInfo;
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.traceInfo.pixelIdxList=traceInfo{iIn}.afterBreakTrace.pixelIdxList;
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.leafInfo=[];
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.state=[];
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.isNode=true;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.nodeInfo=[breakNodeInfo,1];
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.rootInfo=[];
            end
        end
    end
end
end
function nodeInfoIn=getInputInfo(bioTree,nodeInfo)
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
        nodeInfo_pre=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
        nodeInfoIn{iIn}.isNode=true;
        nodeInfoIn{iIn}.nodeInfo= nodeInfo_pre;
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
        rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
        nodeInfoIn{iIn}.isNode=false;
        nodeInfoIn{iIn}.rootInfo= rootInfo;
    end
end
end
function bioTree=removeType1Node(bioTree,nodeList)
newList=[];
countNode=0;
iframe=nodeList(1,1);
for iList=1:size(nodeList,1)
    nodeInfo=nodeList(iList,:);
    bioTree{nodeInfo(1)}.node{nodeInfo(2)}=[];
end
for iNode=1:size(bioTree{iframe}.node,2)
    if ~isempty(bioTree{iframe}.node{iNode})
        countNode=countNode+1;
        newList=[newList;[iNode,countNode]];
    end
end

if ~isempty(newList)
    for iList=1:size(newList,1)
        newNode{newList(iList,2)}=bioTree{iframe}.node{newList(iList,1)};
    end
    bioTree{iframe}.node=newNode;
else
    bioTree{iframe}.node=[];
end
if ~isempty(bioTree{iframe}.node)
    for iNode=1:size(bioTree{iframe}.node,2)
        for iIn=1:size(bioTree{iframe}.node{iNode}.In,2)
            if bioTree{iframe}.node{iNode}.In{iIn}.isNode==true
                nodeInfopre=bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo;
                bioTree{nodeInfopre(1)}.node{nodeInfopre(2)}.Out{nodeInfopre(3)}.nodeInfo=[iframe,iNode,iIn];
            end
            if bioTree{iframe}.node{iNode}.In{iIn}.isNode==false
                rootInfopre=bioTree{iframe}.node{iNode}.In{iIn}.rootInfo;
                bioTree{rootInfopre(1)}.root{rootInfopre(2)}.nodeInfo=[iframe,iNode,iIn];
            end
        end
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
                nodeInfonext=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
                bioTree{nodeInfonext(1)}.node{nodeInfonext(2)}.In{nodeInfonext(3)}.nodeInfo=[iframe,iNode,iOut];
            end
            if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==false
                leafInfonext=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
                bioTree{leafInfonext(1)}.leavies{leafInfonext(2)}.nodeInfo=[iframe,iNode,iOut];
            end
        end
    end
end
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
pixelIdxList2=xresult2+(yresult2-1)*(xMax-xMin+1);
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
xyMin=[xMin,yMin];
BWImageGain(2:end-1,2:end-1)=BWImage;
% BWImageGain=imfill(BWImageGain,'holes');
end
function pixelIdxListIn=getInputMask(bioTree,nodeInfo)
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
        nodeInfo_pre=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
        pixelIdxListIn{iIn}=bioTree{nodeInfo_pre(1)}.node{nodeInfo_pre(2)}.Out{nodeInfo_pre(3)}.traceInfo.pixelIdxList{end};
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
        rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
        if rootInfo(1)==nodeInfo(1)
            pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.rootPixelDetail;
        end
        if  rootInfo(1)<nodeInfo(1)
            pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{end};
        end
    end
end
end

%%% type2 node reduction
function bioTree=type2NodeReduction(bioTree,strength)
nodeCount=0;
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        nodeList=type2NodeinFrame(bioTree{iframe},iframe,bioTree);
        if ~isempty(nodeList)
            for iList=1:size(nodeList,1)
                nodeInfoIn=getInputInfo(bioTree,nodeList(iList,:));
                [traceInfo,canDivideorNot]=fullNodeType2Tracking(bioTree,nodeList(iList,:),strength);
                nodeList(iList,4)=canDivideorNot;
                if canDivideorNot==1
                    bioTree=type2NodeLinker(bioTree,nodeInfoIn,traceInfo);
                end
                nodeCount=nodeCount+1;
            end
            canDivideNode=nodeList(:,4);
            nodeList(canDivideNode==0,:)=[];
            if ~isempty(nodeList)
                bioTree=removeType2Node(bioTree,nodeList);
            end
        end
    end
end
end
function nodeList=type2NodeinFrame(bioTreeFrame,iframe,bioTree)
nodeList=[];
if min(size(bioTreeFrame.node))~=0
    for iNode=1:size(bioTreeFrame.node,2)
        if size(bioTreeFrame.node{iNode}.In,2)==size(bioTreeFrame.node{iNode}.Out,2) && size(bioTreeFrame.node{iNode}.In,2)<=8
            nodeInfo=[iframe,iNode];
            usefulNode=0;
            pixelIdxListIn=getInputMask(bioTree,nodeInfo);
            for iIn=1:size(pixelIdxListIn,2)
                [~,BWImage]=idx2Xy(pixelIdxListIn{iIn},bioTree{1}.imageSize);
                CC=bwconncomp(BWImage);
                regionNum=CC.NumObjects;
                if regionNum>=2
                    usefulNode=1;
                    break
                end
            end
            if usefulNode==0
                nodeList=[nodeList;[iframe,iNode]];
            end
        end
    end
end
end
function bioTree=type2NodeLinker(bioTree,nodeInfoIn,traceInfo)
for iIn=1:size(nodeInfoIn,2)
    if nodeInfoIn{iIn}.isNode==false
        rootInfo=nodeInfoIn{iIn}.rootInfo;
        if traceInfo{iIn}.is2Node==false;
            leafInfo=traceInfo{iIn}.leafInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node=false;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[];
            %             if leafInfo(1)~=rootInfo(1)
            bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
            %             end
            bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=false;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=rootInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=traceInfo{iIn}.traceInfo.pixelIdxList{end};
        end
        if traceInfo{iIn}.is2Node==true;
            nodeInfoT=traceInfo{iIn}.nodeInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=nodeInfoT;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
            bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=false;
            bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=[];
            bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.rootInfo=rootInfo;
        end
    end
    if nodeInfoIn{iIn}.isNode==true
        nodeInfo=nodeInfoIn{iIn}.nodeInfo;
        if traceInfo{iIn}.is2Node==false;
            leafInfo=traceInfo{iIn}.leafInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node=false;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[];
            %             if leafInfo(1)~=nodeInfo(1)
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
            %             end
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=true;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=nodeInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=traceInfo{iIn}.traceInfo.pixelIdxList{end};
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo= leafInfo;
        end
        if traceInfo{iIn}.is2Node==true;
            nodeInfoT=traceInfo{iIn}.nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo= nodeInfoT;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
            bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=true;
            bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=nodeInfo;
        end
    end
end
end
function bioTree=removeType2Node(bioTree,nodeList)
newList=[];
countNode=0;
iframe=nodeList(1,1);
for iList=1:size(nodeList,1)
    nodeInfo=nodeList(iList,:);
    bioTree{nodeInfo(1)}.node{nodeInfo(2)}=[];
end
for iNode=1:size(bioTree{iframe}.node,2)
    if ~isempty(bioTree{iframe}.node{iNode})
        countNode=countNode+1;
        newList=[newList;[iNode,countNode]];
    end
end

if ~isempty(newList)
    for iList=1:size(newList,1)
        newNode{newList(iList,2)}=bioTree{iframe}.node{newList(iList,1)};
    end
    bioTree{iframe}.node=newNode;
else
    bioTree{iframe}.node=[];
end
if ~isempty(bioTree{iframe}.node)
    for iNode=1:size(bioTree{iframe}.node,2)
        for iIn=1:size(bioTree{iframe}.node{iNode}.In,2)
            if bioTree{iframe}.node{iNode}.In{iIn}.isNode==true
                nodeInfopre=bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo;
                bioTree{nodeInfopre(1)}.node{nodeInfopre(2)}.Out{nodeInfopre(3)}.nodeInfo=[iframe,iNode,iIn];
            end
            if bioTree{iframe}.node{iNode}.In{iIn}.isNode==false
                rootInfopre=bioTree{iframe}.node{iNode}.In{iIn}.rootInfo;
                bioTree{rootInfopre(1)}.root{rootInfopre(2)}.nodeInfo=[iframe,iNode,iIn];
            end
        end
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
                nodeInfonext=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
                bioTree{nodeInfonext(1)}.node{nodeInfonext(2)}.In{nodeInfonext(3)}.nodeInfo=[iframe,iNode,iOut];
            end
            if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==false
                leafInfonext=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
                bioTree{leafInfonext(1)}.leavies{leafInfonext(2)}.nodeInfo=[iframe,iNode,iOut];
            end
        end
    end
end
end

%%% type3 node reduction
function bioTree=type3NodeReduction(bioTree,threShold)
nodeCount=0;
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        nodeList=type3NodeinFrame(bioTree{iframe},iframe);
        if ~isempty(nodeList)
            for iList=1:size(nodeList,1)
                nodeInfoIn=getInputInfo(bioTree,nodeList(iList,:));
                [traceInfo,inNum,outNum,canDivideorNot]=fullNodeType3Tracking(bioTree,nodeList(iList,:),threShold);
                if canDivideorNot==1
                    bioTree=type3NodeLinker(bioTree,nodeInfoIn,traceInfo,inNum,nodeList(iList,:));
                    bioTree=removeType3Node(bioTree,nodeList(iList,:),inNum,outNum);
                end
                nodeCount=nodeCount+1;
                %                 disp(nodeCount);
            end
        end
    end
end
end
function nodeList=type3NodeinFrame(bioTreeFrame,iframe)
nodeList=[];
if min(size(bioTreeFrame.node))~=0
    for iNode=1:size(bioTreeFrame.node,2)
        if (size(bioTreeFrame.node{iNode}.In,2)==3 && size(bioTreeFrame.node{iNode}.Out,2)==2)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==2 && size(bioTreeFrame.node{iNode}.Out,2)>=3)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==4 && size(bioTreeFrame.node{iNode}.Out,2)==3)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==3 && size(bioTreeFrame.node{iNode}.Out,2)==4)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==3 && size(bioTreeFrame.node{iNode}.Out,2)==3)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==4 && size(bioTreeFrame.node{iNode}.Out,2)==5)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==5 && size(bioTreeFrame.node{iNode}.Out,2)==4)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==4 && size(bioTreeFrame.node{iNode}.Out,2)==2)
            nodeList=[nodeList;[iframe,iNode]];
        end
    end
end
end
function bioTree=type3NodeLinker(bioTree,nodeInfoIn,traceInfo,inNum,nodeInfo)
if ~(numel(inNum)==1 && inNum==0)
    for iIn=1:numel(inNum)
        if nodeInfoIn{inNum(iIn)}.isNode==false
            rootInfo=nodeInfoIn{inNum(iIn)}.rootInfo;
            if traceInfo{iIn}.is2Node==false;
                leafInfo=traceInfo{iIn}.leafInfo;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node=false;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[];
                %                 if leafInfo(1)~=rootInfo(1)
                bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                %                 end
                bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=leafInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=false;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=rootInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=traceInfo{iIn}.traceInfo.pixelIdxList{end};
            end
            if traceInfo{iIn}.is2Node==true;
                nodeInfoT=traceInfo{iIn}.nodeInfo;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=nodeInfoT;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=false;
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=[];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.rootInfo=rootInfo;
            end
        end
        if nodeInfoIn{inNum(iIn)}.isNode==true
            nodeInfo=nodeInfoIn{inNum(iIn)}.nodeInfo;
            if traceInfo{iIn}.is2Node==false;
                leafInfo=traceInfo{iIn}.leafInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node=false;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[];
                %                 if leafInfo(1)~=nodeInfo(1)
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                %                 end
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=leafInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=true;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=nodeInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=traceInfo{iIn}.traceInfo.pixelIdxList{end};
            end
            if traceInfo{iIn}.is2Node==true;
                nodeInfoT=traceInfo{iIn}.nodeInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo= nodeInfoT;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=true;
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=nodeInfo;
            end
        end
    end
end
if numel(inNum)==1 && inNum==0
    rootNum=size(bioTree{1,nodeInfo(1)}.root,2);
    bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.is2Node=traceInfo{1,1}.is2Node;
    if bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.is2Node==0
        leafInfo=traceInfo{1,1}.leafInfo;
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.leafInfo=leafInfo;
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.nodeInfo=[];
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.rootPixelDetail=traceInfo{1,1}.traceInfo.pixelIdxList{1};
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.traceInfo.pixelIdxList=traceInfo{1,1}.traceInfo.pixelIdxList;
        bioTree{1,leafInfo(1)}.leavies{1,leafInfo(2)}.rootInfo=[nodeInfo(1),rootNum+1];
        bioTree{1,leafInfo(1)}.leavies{1,leafInfo(2)}.nodeInfo=[];
        bioTree{1,leafInfo(1)}.leavies{1,leafInfo(2)}.is2Node=0;
    end
    if bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.is2Node==1
        newNodeInfo=traceInfo{1,1}.nodeInfo;
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.nodeInfo=newNodeInfo;
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.leafInfo=[];
        bioTree{1,newNodeInfo(1)}.node{1,newNodeInfo(2)}.In{1,newNodeInfo(3)}.isNode=0;
        bioTree{1,newNodeInfo(1)}.node{1,newNodeInfo(2)}.In{1,newNodeInfo(3)}.nodeInfo=[];
        bioTree{1,newNodeInfo(1)}.node{1,newNodeInfo(2)}.In{1,newNodeInfo(3)}.rootInfo=[nodeInfo(1),rootNum+1];
    end
    bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.rootPixelDetail=traceInfo{1,1}.traceInfo.pixelIdxList{1};
    bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.traceInfo.pixelIdxList=traceInfo{1,1}.traceInfo.pixelIdxList;
end
end
function bioTree=removeType3Node(bioTree,nodeInfo,inNum,outNum)
if ~(numel(inNum)==1 && inNum==0)
    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In(inNum)=[];
    for i=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,i}.isNode==0
            rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,i}.rootInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo(3)=i;
        end
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,i}.isNode==1
            NewnodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,i}.nodeInfo;
            bioTree{NewnodeInfo(1)}.node{NewnodeInfo(2)}.Out{1,NewnodeInfo(3)}.nodeInfo(3)=i;
        end
    end
end
bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out(outNum)=[];
for i=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,i}.is2Node==0
        leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,i}.leafInfo;
        bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo(3)=i;
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,i}.is2Node==1
        NewnodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,i}.nodeInfo;
        bioTree{NewnodeInfo(1)}.node{NewnodeInfo(2)}.In{1,NewnodeInfo(3)}.nodeInfo(3)=i;
    end
end
end

%%% type4 node reduction
function bioTree=type4NodeReduction(bioTree)
nodeCount=0;
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        nodeList=type4NodeinFrame(bioTree{iframe},iframe);
        if ~isempty(nodeList)
            for iList=1:size(nodeList,1)
                %nodeInfoIn=getInputInfo(bioTree,nodeList(iList,:));
                nodeInfo=nodeList(iList,:);
                pixelIdxListIn=getInputMask(bioTree,nodeInfo);
                iInProps=0;
                for iIn=1:size(pixelIdxListIn,2)
                    [~,BWImage]=idx2Xy(pixelIdxListIn{iIn},bioTree{1}.imageSize);
                    CC=bwconncomp(BWImage);
                    regionNum=CC.NumObjects;
                    if regionNum>=2
                        iInProps=1;
                    end
                end
                if iInProps==1
                    break
                end
                imageSize=bioTree{1}.imageSize;
                traceInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.pixelIdxList;
                pictureOut=false(imageSize);
                pictureOut(traceInfo{1})=true;
                cc=bwconncomp(pictureOut);
                if cc.NumObjects==2
                    outPixelIdxList=cc.PixelIdxList;
                    if numel(traceInfo)>1
                        nodeNum=size(bioTree{nodeInfo(1)+1}.node,2);
                        bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1};
                        bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.traceInfo.pixelIdxList(1)=[];
                        if bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.is2Node==1
                            newNodeInfo=bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.nodeInfo;
                            bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{newNodeInfo(3)}.nodeInfo=[nodeInfo(1)+1,nodeNum+1,1];
                        end
                        if bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.is2Node==0
                            newLeafInfo=bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.leafInfo;
                            bioTree{newLeafInfo(1)}.leavies{newLeafInfo(2)}.nodeInfo=[nodeInfo(1)+1,nodeNum+1,1];
                        end
                        for iOut=1:2
                            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.is2Node=1;
                            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.leafInfo=[];
                            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.nodeInfo=[nodeInfo(1)+1,nodeNum+1,iOut];
                            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.traceInfo.pixelIdxList=[];
                            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.traceInfo.pixelIdxList{1}=outPixelIdxList{iOut};
                        end
                        for iIn=1:2
                            bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.In{1,iIn}.isNode=1;
                            bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.In{1,iIn}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                        end
                    end
                    if numel(traceInfo)==1
                        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.is2Node==1
                            newNodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.nodeInfo;
                            inNum=size(bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In,2);
                            for iIn=1:2
                                if iIn==2
                                    bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,iIn+inNum-1}.isNode=1;
                                    bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,iIn+inNum-1}.rootInfo=[];
                                    bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,iIn+inNum-1}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                                else
                                    bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,newNodeInfo(3)}.isNode=1;
                                    bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,newNodeInfo(3)}.rootInfo=[];
                                    bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,newNodeInfo(3)}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                                end
                            end
                            for iOut=1:2
                                if iOut==2
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.is2Node=1;
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.leafInfo=[];
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.nodeInfo=[newNodeInfo(1),newNodeInfo(2),iOut+inNum-1];
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.traceInfo.pixelIdxList{1}=outPixelIdxList{iOut};
                                else
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.is2Node=1;
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.leafInfo=[];
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.nodeInfo=newNodeInfo;
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.traceInfo.pixelIdxList{1}=outPixelIdxList{iOut};
                                end
                            end
                        end
                    end
                end
                nodeCount=nodeCount+1;
                %             disp(nodeCount)
            end
        end
    end
end
fprintf('\n')
end
function nodeList=type4NodeinFrame(bioTreeFrame,iframe)
nodeList=[];
for iNode=1:size(bioTreeFrame.node,2)
    if size(bioTreeFrame.node{iNode}.In,2)==1 && size(bioTreeFrame.node{iNode}.Out,2)==1
        nodeList=[nodeList;[iframe,iNode]];
    end
end
end
