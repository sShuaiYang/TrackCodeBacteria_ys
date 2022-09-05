function agroseBottomTracking_cBFImages()
dirAllFile='E:\2019-12-23 PAO1_IP32_100x ys';
allNameList=dir(dirAllFile);
%%
for iField=31:numel(allNameList)-2
    if ~(strcmp(allNameList(iField+2).name(1:5),'field') ) && strcmp(allNameList(iField+2).name(end-2:end),'mat')
        continue
    end
    clc
    fieldNum=str2num(allNameList(iField+2).name(end-3:end));
    dirFile=[dirAllFile,'\',allNameList(iField+2).name];
%     dirTree=[dirFile,'\bioTreeResult'];
    
    dirImage=[dirFile,'\Tracking'];
    nameList=dir(dirImage);
    disp(['loading',32,allNameList(iField+2).name]);
    
    %% 生成bioTree    
    tic;
    maskImage=zeros(2048,2048,numel(nameList)-3);    
    parfor i=1:numel(nameList)-3
%         if strcmp(nameList(i+3).name(1:5),'image')
            temp=load([dirImage,'\',nameList(i+3).name]);
            maskImage(:,:,i)=temp.imageTracking;
            %                 maskImage=imerode(maskImage,ones(4));
%             maskImage(:,:,i)=imerode(maskImage(:,:,i),ones(2));%ys 明场相关图像改为ones(2),否则腐蚀太过
            maskImage(:,:,i)=bwmorph(maskImage(:,:,i),'open');%ys 先腐蚀后膨胀，         
%             
%         end
    end    
    maskImages=zeros(size(maskImage,1),size(maskImage,2),size(maskImage,3)*2);
    maskImages(:,:,1:2:end)=maskImage; %jzy用途 1，1，2，2,...构建明确的连通区域
    maskImages(:,:,2:2:end)=maskImage;
    maskImages=logical(maskImages);
    toc;
    
    %         bacNum=getBasicFigure(maskImages,dirFile,mip,fieldNum);
    bacNum=getBasicFigureBF(maskImages,dirFile,fieldNum);
    
%     % 根据空间距离判定到底有几块细菌
%     image1=maskImages(:,:,1);
%     cc=regionprops(image1,'Centroid');
%     if numel(cc)~=0
%         for i=1:numel(cc)
%             info(i,:)=cc(i).Centroid;
%         end
%         if size(info,1)~=1
%             z=linkage(info);
%             %             c=cluster(z,'cutoff',100,'criterion','distance');
%             c=cluster(z,'cutoff',100,'criterion','distance');
%             bacNum=bacNum/max(c);
%         end
%     end
%     cc=bwconncomp(bacNum<=200);
%     indexMask=cc.PixelIdxList{1}(end);
%     maskImages=maskImages(:,:,1:2*indexMask);
    
    disp('treeProcessing')
 
    
%     [maskImages,bestPositionAccumulation]=fluoImageCorrection(maskImages);%xyshift correction    
    dirSave=strcat(dirFile,'\maskImages');   %put the folder of the save results
    mkdir(dirSave)
    tempMasks=maskImages(:,:,1:2:end);
    % save maskImages
    n=0;
    j=0;
    for i=1:size(tempMasks,3)
        n=n+1;
        maskImageStack(:,:,n)=tempMasks(:,:,i);
        if n==200||(j==floor(size(tempMasks,3)/200)&&n==mod(size(tempMasks,3),200))
            j=j+1;
            save([dirSave,'\',num2str(j,'%05.f')],'maskImageStack');
            n=0;
            clear maskImageStack
        end
    end
    
    [~,bioTree]=treeTrackingJob(maskImages,maskImages,1,size(maskImages,3),0);
    
    bioTree{1}.imageSize=[size(maskImages,1),size(maskImages,2)];
%     bioTree{1}.bestPositionAccumulation=bestPositionAccumulation;
    bioTree=bioTreeSizeReduction(bioTree,0);
    bioTree=autoTidyBioTree(bioTree);
%     save([dirFile,'\bioTree_noMearsure.mat'],'bioTree','-v7.3')
    try
        bioTree=agaroseBottomAutoSeg(bioTree);
    catch err
    end
    n=0;
    while n==0
        try
            [bioTree,~,~,~]=myBiograph_new2(bioTree);
            n=1;
        catch err
            bioTree=bioTreeCut(bioTree,numel(bioTree)-2);
        end
    end
    try
        bioTree=type1NodeReduction(bioTree,3,40);
        bioTree=type2NodeReduction(bioTree,2);
        bioTree=type3NodeReduction(bioTree,30);
        bioTree=type4NodeReduction(bioTree);
    catch err
    end
    %%  bioTreeMeasure
    bioTree{1}.fieldNum=fieldNum;
    bioTree=bioTreeDigHoleAndMeasure(bioTree,1,bioTree{1}.imageSize(1),bioTree{1}.imageSize(2));
%     save([dirFile,'\bioTree.mat'],'bioTree','-v7.3')
    % detect cutNum for each branch
    bioTree=detectCutNumForBranch(bioTree);
    bioTree=bioTreeTwoPointAging(bioTree);
    testResult=[dirFile,'\testResult'];
    mkdir(testResult)
    load([dirFile,'\','bacInfo&tree']);
    
    branchList=bioTree{1}.branchList;
    if branchList(1,1)==0
        allData.maskResult=[];
        allData.maskResultTag=[];
        allData.maskResultTagValue=[];
        allData.trackingResult=[];
        allData.trackingResultTag=[];
        allData.trackingResultTagValue=[];
        allData.treeAll=[];
        allData.treeSize=[];
        save([dirFile,'\allData.mat'],'allData')
    else
        [ bioTree ] = bioTreeAllInfoGet_cBF( bioTree,dirFile);
        tree=smallTreeStructureSeperate_cBF(bioTree);
        tree=treeSizeAndGenerationMarker(tree,bioTree{1}.bioTreeTimer(end));
   
        
        %     createAllData(dirFile);
        for iTree=1:numel(tree)
            if tree(iTree).leafNum>=2% original 10
                plotLinkTree(tree(iTree).linkMatrix,tree(iTree).timer,1:tree(iTree).leafNum,[1,0,0],'normal',tree(iTree).cutNum);
                saveas(gcf,[testResult,'\',num2str(iTree,'%03.f'),'.tif'])
                saveas(gcf,[testResult,'\',num2str(iTree,'%03.f'),'.fig'])
                close all
            end
        end
    end
    
    [growthData,~]=generationGrowthRateCal_agrose_new(bioTree,dirFile);
%     growthDataAll{iField}=growthData;
    save([dirFile,'\','bacInfo&tree'],'bioTree','tree','bacInfo','growthData','-v7.3');
end
% save([dirAllFile,'\growthDataAll.mat'],'growthDataAll');
end

%% 作出基本的细菌数目变化图
function bacNum=getBasicFigureBF(maskImages,dirFile,fieldNum)

if isempty(maskImages)
    bacInfo{1}.majorLength=1;
    bacInfo{1}.minorLength=1;
    bacInfo{1}.time=1;
    if ~isempty(bacInfo)
        bacInfo{1}.tag=fieldNum;
    end
    save([dirFile,'\bacInfo&tree'],'bacInfo');
    bacNum=[];
    return
end
load([dirFile,'\Tracking\frameInfo.mat']);
timeBegin=frameInfo(1,1:6);
for iTime=1:size(frameInfo,1)
    timeAll(iTime)=etime(frameInfo(iTime,1:6),timeBegin)/60;
end
maskImages=maskImages(:,:,2:2:end);
dirBasicFigureFile=[dirFile,'\BasicFigure'];
mkdir(dirBasicFigureFile)

% plot bacNum vs time
for iImage=1:size(maskImages,3)
    image=maskImages(:,:,iImage);
    cc=bwconncomp(image);
    bacNum(iImage)=cc.NumObjects;
    cc=regionprops(image,'MajorAxisLength','MinorAxisLength');
    majorLength=[];
    minorLength=[];
    for iSmallBac=1:numel(cc)
        majorLength=[majorLength;cc(iSmallBac).MajorAxisLength];
        minorLength=[minorLength;cc(iSmallBac).MinorAxisLength];
    end
    bacInfo{iImage}.majorLength=majorLength;
    bacInfo{iImage}.minorLength=minorLength;
    bacInfo{iImage}.time=timeAll(iImage);
end
if ~isempty(maskImages)
    plot(timeAll/60,bacNum)
    saveas(gcf,[dirBasicFigureFile,'\bacNum','.tif'])
    saveas(gcf,[dirBasicFigureFile,'\bacNum','.fig'])
    close all
end

save([dirFile,'\bacInfo&tree'],'bacInfo');
end
%%
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
%% detect cutNum for eachBranch
function bioTree=detectCutNumForBranch(bioTree)
n=0;
while n==0
    try
        [bioTree,branchList,~,~]=myBiograph_new2(bioTree);
        branchList(:,5)=0;
        n=1;
    catch err
        bioTree=bioTreeCut(bioTree,numel(bioTree)-2);
    end
end
wrongBranch=[];
for iCoreBranch=1:size(branchList(branchList(:,3)==1,:),1)
    iBranchInfo=branchList(iCoreBranch,:);
    allNode=bioTree{iBranchInfo(1)}.node{iBranchInfo(2)}.allNode;
    [~,order]=sort(allNode(:,1));
    allNode=allNode(order,:);
    branchList(iCoreBranch,5)=numel(bioTree);
    for iNode=1:size(allNode,1)
        nodeInfo=allNode(iNode,:);
        isDivision=isDivsionNode(bioTree,nodeInfo);
        if isDivision==0;
            branchList(iCoreBranch,5)=nodeInfo(1)-1;
            if branchList(iCoreBranch,5)<branchList(iCoreBranch,1)
                wrongBranch=[wrongBranch;iCoreBranch];
            end
            break
        end
    end
end
branchList(wrongBranch,:)=[];
emptyNum=[];
for iBranch=1:size(branchList,1)
    if branchList(iBranch,5)==0
        if branchList(iBranch,1)>numel(bioTree)
            emptyNum=[emptyNum;iBranch];
        else
            branchList(iBranch,5)=numel(bioTree);
        end
    end
end
branchList(emptyNum,:)=[];
bioTree{1}.branchList=branchList;
end
function isDivision=isDivsionNode(bioTree,nodeInfo)
isDivision=0;
if size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)==1 && size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)==2
    %     frameOut=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out;
    %     bacteriaSize=[frameOut{1}.traceInfo.measurment{1}.FilledArea,frameOut{2}.traceInfo.measurment{1}.FilledArea];
    %     if ~(max(bacteriaSize)-min(bacteriaSize)>150 || max(bacteriaSize)>min(bacteriaSize)*2 )
    isDivision=1;
    %     end
end
end
%%
function tree=treeSizeAndGenerationMarker(tree,calTime)
for iTree=1:numel(tree)
    try
        smallTree=tree(iTree);
        linkInfo=smallTree.linkInfo;
        tree(iTree).linkInfo(1).linkGeneration=1;
        linkGeneration{1}=1;
        linkGeneration{2}=1;
        for iLink=2:numel(linkInfo)
            targetInfo=smallTree.linkColumn(smallTree.linkRow==smallTree.linkRow(iLink));
            linkTagetNum=numel(targetInfo);
            index=find(targetInfo==smallTree.linkColumn(iLink));
            linkGeneration{smallTree.linkColumn(iLink)}=[linkGeneration{smallTree.linkRow(iLink)};index];
            linkInfo(iLink).linkGeneration=linkGeneration{smallTree.linkColumn(iLink)};
            tree(iTree).linkInfo(iLink).linkGeneration=numel(linkInfo(iLink).linkGeneration);
        end
    catch err
        linkInfo=smallTree.linkInfo;
        for iLink=1:numel(linkInfo)
            tree(iTree).linkInfo(iLink).linkGeneration=NaN;
        end
    end
    
end
end
%% xyshift correction
function [imageStack,bestPositionAccumulation] = fluoImageCorrection(imageStack)
% 针对荧光图像的时间序列的xy矫正
% 输入变量为二值图像，核心矫正函数已经优化过
imageSize=size(imageStack);
if min(imageSize(1:2))>=1500
    imageStackNew=imageStack(1024-500:1024+500,1024-500:1024+500,:);
%     imageStackNew=imageStack(600:1100,600:1100,:);
else
    imageStackNew=imageStack;
end
imageStackNew1(:,:,2:size(imageStackNew,3)+1)=imageStackNew;
parfor i=2:size(imageStack,3)
    bestPosition(i,:)=caculateCrossCorrelationForImage(imageStackNew1(:,:,i),imageStackNew(:,:,i),15);
end
imageStack=imageCorrectionWithBestPosition(imageStack,bestPosition);
bestPositionAccumulation=zeros(size(bestPosition));
for i=2:size(bestPosition,1)
    bestPositionAccumulation(i,:)=bestPosition(i,:)+bestPositionAccumulation(i-1,:);
end
% gfpImage=imageCorrectionWithBestPosition(gfpImage,bestPosition);
% rfpImage=imageCorrectionWithBestPosition(rfpImage,bestPosition);
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
%%
