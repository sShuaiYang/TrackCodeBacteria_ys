function bioTree=velocityTracking(bioTree)
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        nodeInfo=[];
        isNode=false;
        for iroot=1:size(bioTree{iframe}.root,2)
            rootInfo=[iframe,iroot];
           [bioTree, dataList]=getData(bioTree,rootInfo,nodeInfo,isNode);
            if ~isempty(dataList)
                [dataListDenoise,dataVelocity] =getVelocity(dataList);
            else
                dataVelocity=[];
            end
            if ~isempty(dataVelocity)
                bioTree=setData(bioTree,rootInfo,nodeInfo,isNode,dataListDenoise,dataVelocity);
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        rootInfo=[];
        isNode=true;
        for iNode=1:size(bioTree{iframe}.node,2)
            for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                nodeInfo=[iframe,iNode,iOut];
                [bioTree,dataList]=getData(bioTree,rootInfo,nodeInfo,isNode);
                if ~isempty(dataList)
                    [dataListDenoise,dataVelocity] =getVelocity(dataList);
                else
                    dataVelocity=[];
                end
                if ~isempty(dataVelocity)
                    bioTree=setData(bioTree,rootInfo,nodeInfo,isNode,dataListDenoise,dataVelocity);
                end
            end
        end
    end
end
end
function [bioTree, dataList]=getData(bioTree,rootInfo,nodeInfo,isNode) %calulate vector Orirntation and return the data List in each root or a node% dataList=[p1xPosition',p1yPosition',p2xPosition',p2yPosition',centroidx',centroidy',majorAxisLength'];
dataList=[];
p1xPosition=[];p1yPosition=[];
p2xPosition=[];p2yPosition=[];
centroidx=[];centroidy=[];
majorAxisLength=[];

if isNode==false
    if ~isempty(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment)
        for iTrace=1:size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment,2)
            p1_x=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).p1Position(1);
            p1_y=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).p1Position(2);
            p2_x=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).p2Position(1);
            p2_y=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).p2Position(2);
            p1xPosition=[ p1xPosition,p1_x];
            p1yPosition=[ p1yPosition,p1_y];
            p2xPosition=[ p2xPosition,p2_x];
            p2yPosition=[ p2yPosition,p2_y];
            centroidx=[ centroidx,bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).Centroid(1)];
            centroidy=[ centroidy,bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).Centroid(2)];
            majorAxisLength=[majorAxisLength,bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).MajorAxisLength];
        end
        dataList=[p1xPosition',p1yPosition',p2xPosition',p2yPosition',centroidx',centroidy',majorAxisLength'];
        return;
    end
end
if isNode==true
    if ~isempty(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment)
        for iTrace=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment,2)
            p1_x=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).p1Position(1);
            p1_y=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).p1Position(2);
            p2_x=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).p2Position(1);
            p2_y=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).p2Position(2);
            p1xPosition=[ p1xPosition,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).p1Position(1)];
            p1yPosition=[ p1yPosition,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).p1Position(2)];
            p2xPosition=[ p2xPosition,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).p2Position(1)];
            p2yPosition=[ p2yPosition,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).p2Position(2)];
            centroidx=[ centroidx,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).Centroid(1)];
            centroidy=[ centroidy,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).Centroid(2)];
            majorAxisLength=[majorAxisLength,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).MajorAxisLength];
        end
        dataList=[p1xPosition',p1yPosition',p2xPosition',p2yPosition',centroidx',centroidy',majorAxisLength'];
        return;
    end
end
end
function [dataListDenoise,dataVelocity] =getVelocity(dataList)
dataListDenoise=zeros(size(dataList,1),size(dataList,2));
level=4;
for iData=1:size(dataList,2)
    dataListDenoise(:,iData)=wden(dataList(:,iData),'sqtwolog','s','sln',level,'sym8');
end
if size(dataListDenoise,1)>1
    dataVelocityTemp=zeros(size(dataList,1)-1,size(dataList,2));
    for iData=1:size(dataList,2)
        dataVelocityTemp(:,iData)=diff(dataListDenoise(:,iData));
    end
    p1VelocityAmp=sqrt(dataVelocityTemp(:,1).^2+dataVelocityTemp(:,2).^2);
    p2VelocityAmp=sqrt(dataVelocityTemp(:,3).^2+dataVelocityTemp(:,4).^2);
    centroidVelocityAmp=sqrt(dataVelocityTemp(:,5).^2+dataVelocityTemp(:,6).^2);
    angleV1Pre=angleBetweenTwoVector(dataVelocityTemp(:,1:2),[dataListDenoise(1:end-1,1)-dataListDenoise(1:end-1,3),dataListDenoise(1:end-1,2)-dataListDenoise(1:end-1,4)]);
    angleV2Pre=angleBetweenTwoVector(dataVelocityTemp(:,3:4),[dataListDenoise(1:end-1,1)-dataListDenoise(1:end-1,3),dataListDenoise(1:end-1,2)-dataListDenoise(1:end-1,4)]);
    angleV1next=angleBetweenTwoVector(dataVelocityTemp(:,1:2),[dataListDenoise(2:end,1)-dataListDenoise(2:end,3),dataListDenoise(2:end,2)-dataListDenoise(2:end,4)]);
    angleV2next=angleBetweenTwoVector(dataVelocityTemp(:,3:4),[dataListDenoise(2:end,1)-dataListDenoise(2:end,3),dataListDenoise(2:end,2)-dataListDenoise(2:end,4)]);
    angleVelocity=angleBetweenTwoVector([dataListDenoise(2:end,1)-dataListDenoise(2:end,3),dataListDenoise(2:end,2)-dataListDenoise(2:end,4)],[dataListDenoise(1:end-1,1)-dataListDenoise(1:end-1,3),dataListDenoise(1:end-1,2)-dataListDenoise(1:end-1,4)]);
    dataVelocity(:,1:2)=dataVelocityTemp(:,1:2);
    dataVelocity(:,3)=p1VelocityAmp;
    dataVelocity(:,4)= angleV1Pre;
    dataVelocity(:,5)= angleV1next;
    dataVelocity(:,6:7)=dataVelocityTemp(:,3:4);
    dataVelocity(:,8)=p2VelocityAmp;
    dataVelocity(:,9)= angleV2Pre;
    dataVelocity(:,10)= angleV2next;
    dataVelocity(:,11:12)=dataVelocityTemp(:,5:6);
    dataVelocity(:,13)=centroidVelocityAmp;
    dataVelocity(:,14)=angleVelocity;
    dataVelocity(:,15)=dataVelocityTemp(:,7);
else
    dataVelocity=[];
end
end

function angle=angleBetweenTwoVector(vector1,vector2) %return the angle ([-180,180]) betwwen the vector1 and vector2

angle1=atan(vector1(:,2)./abs(vector1(:,1))).*(180/pi);
angle1(vector1(:,1)<0&vector1(:,2)>=0)= 180-angle1(vector1(:,1)<0&vector1(:,2)>=0);
angle1(vector1(:,1)<0&vector1(:,2)<0)= -180-angle1(vector1(:,1)<0&vector1(:,2)<0);

angle2=atan(vector2(:,2)./abs(vector2(:,1))).*(180/pi);
angle2(vector2(:,1)<0&vector2(:,2)>=0)= 180-angle2(vector2(:,1)<0&vector2(:,2)>=0);
angle2(vector2(:,1)<0&vector2(:,2)<0)= -180-angle2(vector2(:,1)<0&vector2(:,2)<0);

angle=angle2- angle1;
angle(angle>180)=angle(angle>180)-360;
angle(angle<-180)=angle(angle<-180)+360;
end

function  bioTree=setData(bioTree,rootInfo,nodeInfo,isNode,dataListDenoise,dataVelocity)
if isNode==false
    for iTrace=1:size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment,2)-1
        bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).denosieData=dataListDenoise(iTrace,:);
        bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).velocityData=dataVelocity(iTrace,:);
    end
      bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{end}(1).denosieData=dataListDenoise(end,:);
    return;
end
if isNode==true
    for iTrace=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment,2)-1
       bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).denosieData=dataListDenoise(iTrace,:);
       bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).velocityData=dataVelocity(iTrace,:);
    end
     bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{end}(1).denosieData=dataListDenoise(end,:);
    return;
end
end
