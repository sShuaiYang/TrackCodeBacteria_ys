function getBioTreeAndCalculateMiuAndProteinProduceRate1(dirFile,image)
% 针对Jzy & Xag 分析XyScan montage的板底细菌HJD1而设计的批处理程序
% 理论上程序流程分一下步骤
% 1.montage存储的时间序列的图片是分离的，首先进行整合，变成扫描视野个数的时间序列，并进行背景矫正（用罗丹明与钙黄绿素矫正）
% 2.进行图像二值转换，细菌识别
% 3.对图像xy漂移进行矫正，并将结果的红绿图像进行保存
% 4.由maskImage进行制作bioTree，并保存
% 5.从荧光图片中获得对应的荧光信息并计算对应的数值
% 数据结构为  目标文件夹--allResult--imageGFP--imageGFP001,imageGFP002
 
% dirFile='F:\data\2015-11-18 EMCCD HJD1_B test\2015-11-18 jzy test1';
% dirFile=uigetdir();

% 自动将single tif的文件整理好并放入对应的文件夹。
% manualTidyForSingleImage(dirFile);
% return
% montage的整合以及背景矫正
% montageImageGet(dirFile);

% return
dirAllResult=strcat(dirFile,'\allResult');
nameList=dir(dirAllResult);
for iField=1:numel(nameList)-2
    clc
    %%
%     iImage=image(:,:,iField);
    iImage=0;
    disp(iField)
    tic;
    fieldName=strcat(dirAllResult,'\',nameList(iField+2).name);
    cd(fieldName)
    dirMaskImage=strcat(fieldName,'\','maskImage');
    dirBioTreeResult=strcat(fieldName,'\','bioTreeResult');
    dirMiuResult=strcat(fieldName,'\','miu&Presult');
    mkdir(dirMaskImage)
    mkdir(dirBioTreeResult)
    mkdir(dirMiuResult)
    dirGFPImage=strcat(fieldName,'\imageGFP');
    dirRFPImage=strcat(fieldName,'\imageRFP');
% %     这里暂时没有加入批处理，只处理一个gfpImageStack
%     load(strcat(dirGFPImage,'\imageGFP01'));
%     load(strcat(dirRFPImage,'\imageRFP01'));
%     
%     图像的二值转换
%     imageGFP=imageGFP(:,:,1:12);
%     imageRFP=imageRFP(:,:,1:12);   %% add by jzy 11.3
%     maskImage=fluoImageProcessing(imageGFP);
% %     
% %     % xy漂移的矫正   
%     [maskImage,imageGFP,imageRFP]=fluoImageCorrection(maskImage,imageGFP,imageRFP);
% %     
% %     % 将结果储存一下
%      save(strcat(dirMaskImage,'\1'),'maskImage','imageGFP','imageRFP')
%      load(strcat(dirMaskImage,'\1'))
%     
%     % 从maskImage制作bioTree
%     [~,bioTree]=treeTrackingJob([],maskImage,1,size(maskImage,3),0);
%     bioTree{1}.imageSize=size(maskImage(:,:,1));
% %     bioTree=leafRootRefineBioTree(bioTree,0);
% %     bioTree=bioTreeSizeReduction(bioTree,0);
% %     bioTree=autoTidyBioTree(bioTree);
% %     bioTree=type2NodeReduction(bioTree,1);
% %     bioTree=mixMatchReducion(bioTree,5);
% %     bioTree=type3NodeReduction(bioTree,15);
% %     bioTree=type4NodeReduction(bioTree);
% %     bioTree=bioTreeSegmentation(bioTree);
%     bioTree=bioTreeDigHoleAndMeasure(bioTree,1,bioTree{1}.imageSize(1),bioTree{1}.imageSize(2));   % 实心的
%     dirSmallTree=strcat(dirBioTreeResult,'\','1');
%     save(dirSmallTree,'bioTree')
%      
%     load(strcat(dirMaskImage,'\1'))
%     load(dirSmallTree);
%     imageRFP=backGroundSubstract(imageRFP);
%     imageGFP=backGroundSubstract(imageGFP);
%     [bacInfo]=mainInfomationGet(bioTree,imageGFP,imageRFP);
    dirBacInfo=strcat(fieldName,'\','bacInfo');
%     mkdir(dirBacInfo)
%     save(strcat(dirBacInfo,'\','2_选用一些较长的片段用于计算提高正确性'),'bacInfo')
    
    load(strcat(dirBacInfo,'\','2_选用一些较长的片段用于计算提高正确性'));
    [gama,miu,fpProduce,GRratio,gfp,rfp,isDead]=calculateMiuAndProduceRate(bacInfo,iImage);
    save(strcat(dirMiuResult,'\','3_只计算第一个点_new_1点'),'gama','miu','fpProduce','GRratio','gfp','rfp','isDead')
    toc;
end
end

% z之前的主程序
% function getBioTreeAndCalculateMiuAndProteinProduceRate()
% % 针对Jzy & Xag 分析XyScan montage的板底细菌HJD1而设计的批处理程序
% % 理论上程序流程分一下步骤
% % 1.montage存储的时间序列的图片是分离的，首先进行整合，变成扫描视野个数的时间序列，并进行背景矫正（用罗丹明与钙黄绿素矫正）
% % 2.进行图像二值转换，细菌识别
% % 3.对图像xy漂移进行矫正，并将结果的红绿图像进行保存
% % 4.由maskImage进行制作bioTree，并保存
% % 5.从荧光图片中获得对应的荧光信息并计算对应的数值
% % 数据结构为  目标文件夹--allResult--imageGFP--imageGFP001,imageGFP002
%                                   
% dirFile=uigetdir();
% dirMiuResult=strcat(dirFile,'\miu&P_result');
% for iField=55:399
%     disp(iField)
%     dirMaskImage=strcat(dirFile,'\maskImage\',num2str(iField));
%     load(dirMaskImage)
%     dirBioTreeResult=strcat(dirFile,'\bioTreeResult\',num2str(iField));
%     load(dirBioTreeResult)
%     
% %     % 图像的二值转换
% %     maskImage=fluoImageProcessing(imageGFP);
% %     
% %     % xy漂移的矫正
% %     [maskImage,imageGFP,imageRFP]=fluoImageCorrection(maskImage,imageGFP,imageRFP);
% %     
% %     % 将结果储存一下
% %     save(strcat(dirMaskImage,'\',num2str(iField)),'maskImage','imageGFP','imageRFP')
% %     
% %     % 从maskImage制作bioTree
% %     [~,bioTree]=treeTrackingJob([],maskImage,1,size(maskImage,3),0);
% %     bioTree{1}.imageSize=size(maskImage(:,:,1));
% %     bioTree=leafRootRefineBioTree(bioTree,0);
% %     bioTree=bioTreeSizeReduction(bioTree,0);
% % %     bioTree=bioTreeSegmentation(bioTree);
% %     bioTree=bioTreeDigHoleAndMeasure(bioTree,1,bioTree{1}.imageSize(1),bioTree{1}.imageSize(2));   % 实心的
% %     dirSmallTree=strcat(dirBioTreeResult,'\',num2str(iField));
% %     save(dirSmallTree,'bioTree')
% %     
%     imageRFP=backGroundSubstract(imageRFP);
%     imageGFP=backGroundSubstract(imageGFP);
%     [gama,miu,fpProduce]=calculateMiuAndProduceRate(bioTree,imageGFP,imageRFP);
%     save(strcat(dirMiuResult,'\',num2str(iField)),'gama','miu','fpProduce')
% end
% end
%% 图像二值化程序
function maskImage=fluoImageProcessing(imageStack)
parfor iStack=1:size(imageStack,3)
    image=imageStack(:,:,iStack);
    pixelInfo=image(:);
    pixelInfo=sort(pixelInfo);
    warning off
    backGround=mean(pixelInfo(1:end/3));
    image=image-backGround;
    edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
    gaussianFilter1=fspecial('gaussian',[5, 5],20);
    gaussianFilter2=fspecial('gaussian',[3, 3],2);
    image=uint8(double(image)/double(max(max(image)))*200);
    image=imfilter(image,gaussianFilter1); % use Guassian blur filter process
    image=imfilter(image,edgeFilter); %use edgeFilter process
%     image=imfilter(image,gaussianFilter1);
    image=im2bw(image,40/255);
    image=bwareaopen(image,30);
    image=logical(image);
    image=imclearborder(image);
    image=imfill(image,'holes');
    maskImage(:,:,iStack)=image;
%     maskImage(:,:,iStack)=bwmorph(image,'remove');
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

%% backGround substract
function  [gfpImage]=backGroundSubstract(gfpImage)
%首先要对红色荧光和绿色荧光图像去除背景
for i=1:size(gfpImage,3)
    gImage=gfpImage(:,:,i);
    gIndex=gImage(:);
    gIndex=sort(gIndex);
    gBack=mean(gIndex(1:round(end/3)));
    gfpImage(:,:,i)=gfpImage(:,:,i)-gBack;
end
end

%% 提取计算时需要信息的主程序
function [bacInfo]=mainInfomationGet(bioTree,gfpImage,rfpImage)
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
    meanRFP=[];
    for iTrace=1:numel(pixelIdxList);
        pixelInfo=pixelIdxList{iTrace};
%         [xyMin,BWImageGain]=idx2Xy(pixelInfo,[512,512]);
%         pixelInfo=xy2Idx(xyMin,BWImageGain,[512,512]);
        gImage=gfpImage(:,:,stackInfo(iTrace));
        rImage=rfpImage(:,:,stackInfo(iTrace));
        meanGFP(iTrace,1)=mean(double(gImage(pixelInfo)));
        meanRFP(iTrace,1)=mean(double(rImage(pixelInfo)));  
    end
    bacInfo{i}.lengthInfo=pieceInfo{i}.length;
    bacInfo{i}.meanGFP=meanGFP;
    bacInfo{i}.meanRFP=meanRFP;
    bacInfo{i}.time=10*pieceInfo{i}.traceRange;
    bacInfo{i}.time=bacInfo{i}.time';
    bacInfo{i}.pixelInfo=pixelIdxList{end};
end
end

%% 计算各参数及最终数据的主程序
function [gamaAll,miuAll,fpProduceAll,GRratioAll,gfpAll,rfpAll,isDead]=calculateMiuAndProduceRate(bacInfo,image)
% 板底扫描计算细菌的生长速率以及荧光蛋白的生成速率，由原始程序GRratioAnalysis稍加改动而成
% alphaG=0.2840*0.6*1.5712;
% alphaR=0.2878*0.6*3.7625*1.15;   % 1 result
% alphaG=0.4840*0.6;
alphaG=0.3320*0.6;   % 2016-3-8红绿修正
alphaR=1.155*0.6;   % 2 result
% ratio=0.2662;
% E=1-alphaR/alphaG/ratio;
E=0;
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
    return
else
    gamaAll=[];
    miuAll=[];
    fpProduceAll=[];
    GRratioAll=[];
    gfpAll=[];
    rfpAll=[];
    isDead=[];
%     gamaAll=zeros(numel(bacInfo),1);
%     miuAll=zeros(numel(bacInfo),1);
%     fpProduceAll=zeros(numel(bacInfo),1);
end
for i=1:numel(bacInfo)
% %     lengthInfo=bacInfo{i}.lengthInfo;
% % %     if numel(lengthInfo)<12
% % %         continue
% % %     end
% %     meanGFP=bacInfo{i}.meanGFP;
% %     meanRFP=bacInfo{i}.meanRFP;
% %     timeSeries=bacInfo{i}.time;
% %     gama=getGrowthRate(timeSeries,lengthInfo);
% %     if isempty(gama)
% %         continue
% %     end
% %     meanGFP=meanGFP./alphaG;
% %     meanRFP=meanRFP./alphaR;
% %     meanGFP=meanGFP+meanRFP.*E;
% %     meanGFP=smooth(meanGFP);
% %     meanRFP=smooth(meanRFP);
% % %     plot(GRratio)
% % %     meanGFP=sfGFPbleachCorrection(meanGFP);
% % %     meanRFP=rFPbleachCorrection(meanRFP);
% %     GRratio=meanGFP./meanRFP;
% %     if var(GRratio)/mean(GRratio)>0.05
% %        continue
% %     end         
% % %     hold on;plot(meanGFP,'g');plot(meanRFP,'r')
% % %     plot(meanRFP./meanGFP)
% % %     close all
% % %     plot(GRratio)
% %     deltaFP=meanGFP-meanRFP;
% %     deltaT=diff(timeSeries);
% %     deltaT(end+1)=deltaT(end);
% %     deltaGFP=diff(meanGFP);
% %     deltaGFP(end+1)=deltaGFP(end);
% %     deltaRFP=diff(meanRFP);
% %     deltaRFP(end+1)=deltaRFP(end);
% %     
% %     % gama就是由菌长计算出的生长速率，另由红绿比计算一下作为参照，假设mG很大
% %     miu1=(-deltaRFP./deltaT+rMature*deltaFP)./meanRFP;
% % %     miu1=rMature*(meanGFP./meanRFP-1);
% % %     miu1=-(GRratio-1)./(GRratio/gMature-1/rMature);
% %     % gama修正
% %     
% % %     miu2=-((meanGFP-meanRFP)+deltaGFP./deltaT/gMature-deltaRFP./deltaT/rMature)./(meanGFP/gMature-meanRFP/rMature);
% %     
% %     % P为荧光蛋白的生产速率
% %     d_deltaFP=diff(deltaFP);
% %     d_deltaFP(end+1)=d_deltaFP(end);
% %     fpProduce=(gama+rMature).*deltaFP+d_deltaFP./deltaT;
% %     
% % %     plot(fpProduce)
% %     miu1(gama<=0)=[];
% % %     miu(gama<=0)=[];
% %     gama(gama<=0)=[];
% % %     plot(gama);hold on;plot(miu1,'r');
% % %     close 
% %     gamaAll=[gamaAll;mean(gama)];
% %     miuAll=[miuAll;mean(miu1)];
% %     fpProduceAll=[fpProduceAll;mean(fpProduce)];
% %     GRratioAll=[GRratioAll;mean(GRratio)];
% %     gfpAll=[gfpAll;mean(meanGFP)];
% %     rfpAll=[rfpAll;mean(meanRFP)];
% % %     gamaAll(i)=mean(gama);
% % %     miuAll(i)=mean(miu1);
% % %     fpProduceAll(i)=mean(fpProduce);

% changed by jzy 2016.3.8 用第一二张的图片和生长速率来计算
    lengthInfo=bacInfo{i}.lengthInfo;
    timePoint=1;
    if numel(lengthInfo)<timePoint+1
        continue
    end
% %     meanGFP=bacInfo{i}.meanGFP*6510/6510;
% %     meanRFP=bacInfo{i}.meanRFP*11918/9877;     % *11918/9877参考背景11918,实际背景9877
    meanGFP=bacInfo{i}.meanGFP*6510/9000;    %9941
    meanRFP=bacInfo{i}.meanRFP*11918/8812.107;   % 2016.03.20
%     meanRFP=bacInfo{i}.meanRFP*11918/11918; 
    timeSeries=bacInfo{i}.time;
    gama=getGrowthRate(timeSeries,lengthInfo,timePoint);
    if isempty(gama)
        continue
    end
    meanGFP=meanGFP-meanRFP/7.2;
    meanRFP=meanRFP./alphaR;
    meanGFP=meanGFP./alphaG;
%     meanGFP=smooth(meanGFP);
%     meanRFP=smooth(meanRFP);
    if meanGFP(1)<meanRFP(1)
        continue
    end
    GRratio=meanGFP(timePoint)./meanRFP(timePoint);
    deltaFP=meanGFP(timePoint)-meanRFP(timePoint);
    deltaT=diff(timeSeries);
    deltaT=deltaT(1);
    deltaGFP=meanGFP(timePoint+1)-meanGFP(timePoint);
    deltaRFP=meanRFP(timePoint+1)-meanRFP(timePoint);
    
    if deltaGFP<deltaRFP
        continue
    end
    % gama就是由菌长计算出的生长速率，另由红绿比计算一下作为参照，假设mG很大
    miu1=(-deltaRFP./deltaT+rMature*deltaFP)./meanRFP(timePoint);   % 只计算第一个点
    mG=bacInfo{i}.meanGFP*6510/9000;
    mR=bacInfo{i}.meanRFP*11918/8812.107;
    fpProduce=((gama-(mR(timePoint+1)-mR(timePoint))/mR(timePoint)/deltaT)/rMature+1+1/7.2)/(mG(timePoint)/mR(timePoint)*alphaR/alphaG);
    
    if miu1<0
        continue
    end
    
    % P为荧光蛋白的生产速率
%     fpProduce=miu1*meanGFP(timePoint)+deltaGFP/deltaT;
    
%     if any(image(bacInfo{i}.pixelInfo)>=1000)
%         isDead=[isDead;1];
%     else
%         if all(image(bacInfo{i}.pixelInfo)==0)
%             isDead=[isDead;2];
%         else
%             isDead=[isDead;0];
%         end
%     end
    isDead=0;
    gamaAll=[gamaAll;mean(gama)];
    miuAll=[miuAll;mean(miu1)];
    fpProduceAll=[fpProduceAll;mean(fpProduce)];
    GRratioAll=[GRratioAll;mean(GRratio)];
    gfpAll=[gfpAll;mean(meanGFP)];
    rfpAll=[rfpAll;mean(meanRFP)];
%     gamaAll(i)=mean(gama);
%     miuAll(i)=mean(miu1);
%     fpProduceAll(i)=mean(fpProduce);
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
if numel(lengthInfo)>=5 && lengthInfo(end)<lengthInfo(end-1)
    lengthInfo(end)=[];
    timeSeries(end)=[];
end
% gama=prepareCurveData(timeSeries,lengthInfo);
% if gama<=0
%     gama=[];
%     return
% end
% gama=ones(numel(lenT),1)*gama;
% return
% gama(1,1)=(lengthInfo(2)-lengthInfo(1))/(timeSeries(2)-timeSeries(1))/(lengthInfo(1));
i=timePoint;
% while gama<0
%     i=i+1;
%     if i>=numel(lengthInfo)
%         break
%     end
    gama=(lengthInfo(i+1)-lengthInfo(i))/(timeSeries(i+1)-timeSeries(i))/(lengthInfo(i));
% end
if gama<0
    gama=[];
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