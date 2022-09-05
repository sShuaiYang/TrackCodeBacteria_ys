function correctionList=autoGetCorrectionList(bioTree,allDataOri)
dataNum=size(allDataOri,2);
pictureSize=bioTree{1}.imageSize;
correctionList=[];
angleSpace=linspace(-180,180,37);
for iData=1:dataNum
    childCorrectionList=[];
    centroidInfo=allDataOri{iData}.originalData(:,6:7);
    lengthInfo=allDataOri{iData}.originalData(:,9);
    diffLength=abs(diff(lengthInfo));
    diffLength=(diffLength>=10);
    diffLength=[0;diffLength];
    diffLength=1-diffLength;
    centroidInfo(:,1)=centroidInfo(:,1)>40 & centroidInfo(:,1)<=pictureSize(2)-39;
    centroidInfo(:,2)=centroidInfo(:,2)>40 & centroidInfo(:,2)<=pictureSize(1)-39;
    properCentroid=centroidInfo(:,1) & centroidInfo(:,2) & diffLength;
    cc=bwconncomp(properCentroid);
    frameNum=allDataOri{iData}.originalData(:,1);
    for iObject=1:cc.NumObjects
        if numel(cc.PixelIdxList{iObject})>=200
            [recentHead1,recentTail1]=caculateHead(allDataOri{iData}.velocityData(cc.PixelIdxList{iObject},[4,5]),angleSpace);
            %[recentHead2,recentTail2]=caculateHead(allDataOri{iData}.velocityData(cc.PixelIdxList{iObject},[9,10]),angleSpace);
            if recentHead1>recentTail1
                childCorrectionList=[childCorrectionList,iData,frameNum(cc.PixelIdxList{iObject}(1)),frameNum(cc.PixelIdxList{iObject}(end)),0];
            else
                childCorrectionList=[childCorrectionList,iData,frameNum(cc.PixelIdxList{iObject}(1)),frameNum(cc.PixelIdxList{iObject}(end)),1];
            end
        end       
    end
    if ~isempty(childCorrectionList)
        correctionList=catColumnVector(correctionList,childCorrectionList);
    end     
end
     correctionList=findAndDeleteColumn(correctionList);
end
function p1VelocityOrientation=histOrientation(velocityList,angleSpace)
[countP1angle,angle]=getAnglehist(velocityList(:,1),velocityList(:,2),angleSpace);
p1VelocityOrientation.histOrientation=[angle,countP1angle'];
end
function [countP1angle,angle]=getAnglehist(velocityP1,velocityAngle,angleSpace)
countP1angle=[];
angle=[];
for i=2:2:size(angleSpace,2)
    countP1angle=[countP1angle,sum(velocityP1(velocityAngle>=angleSpace(i-1)&velocityAngle<angleSpace(i+1)))];
    angle=[angle,angleSpace(i)];
end
angleRed=angle.*(pi/180);
angle=[angleRed',angle'];
end
function [recentHead,recentTail]=caculateHead(velocityList,angleSpace)
p1VelocityHist=histOrientation(velocityList,angleSpace);
dataList=p1VelocityHist.histOrientation(:,3);
theta=p1VelocityHist.histOrientation(:,2);
recentHead=sum(dataList(theta<90 & theta>-90));
recentTail=sum(dataList(theta<-90 | theta>90));
end
function parentMat=catColumnVector(parentMat,childMat)
lineNum=size(parentMat,1)+1;
parentMat(lineNum,1:size(childMat,2))=childMat;
end
function correctionList=findAndDeleteColumn(correctionList)
columnNum=size(correctionList,2);
deleteColumn=(1:(columnNum-1)/4)*4+1;
correctionList(:,deleteColumn)=[];
end