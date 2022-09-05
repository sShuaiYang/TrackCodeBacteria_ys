function traceResult=bioTreeTraceExplore(bioTree) % this function use to get the velocity of the bacteria
cutOffTime=800;
stepTime=10;
velocityBin=logspace(-3,1.5,80);
timeInte=1:stepTime:cutOffTime;
walkingLengthThreshold=12;
meanDispTime0=zeros(size(timeInte,2),2);
meanDispTime1=zeros(size(timeInte,2),2);
countSum0=zeros(size(velocityBin,2),1);
countSum1=zeros(size(velocityBin,2),1);
walkingTime0=[];
crawlingTime0=[];
walkingTime1=[];
crawlingTime1=[];
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.isCluster==false
                dataList=[];
                for iTrace=1:size(bioTree{iframe}.root{iroot}.traceInfo.measurment,2);
                    dataList=[dataList;[bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(1).Centroid,bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(1).MajorAxisLength]];
                end
                [traceResult,countV,walkingCrawlTime]=traceExploring(dataList,timeInte,velocityBin,walkingLengthThreshold);
                bioTree{iframe}.root{iroot}.traceResult=traceResult;
                bioTree{iframe}.root{iroot}.velocityDistibution=[velocityBin',countV];
                meanDispTime0(:,1)=meanDispTime0(:,1)+traceResult(:,2);
                meanDispTime0(:,2)=meanDispTime0(:,2)+traceResult(:,3);
                countSum0=countSum0+countV;
                walkingTime0=[walkingTime0,walkingCrawlTime(1)];
                crawlingTime0=[crawlingTime0,walkingCrawlTime(2)];
            end
            if bioTree{iframe}.root{iroot}.isCluster==true
                dataList=[];
                for iTrace=1:size(bioTree{iframe}.root{iroot}.traceInfo.measurment,2);
                    dataList=[dataList;[bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(1).Centroid,bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(1).MajorAxisLength]];
                end
                [traceResult,countV,walkingCrawlTime]=traceExploring(dataList,timeInte,velocityBin,walkingLengthThreshold);
                bioTree{iframe}.root{iroot}.traceResult=traceResult;
                bioTree{iframe}.root{iroot}.velocityDistibution=[velocityBin',countV];
                meanDispTime1(:,1)=meanDispTime1(:,1)+traceResult(:,2);
                meanDispTime1(:,2)=meanDispTime1(:,2)+traceResult(:,3);
                countSum1=countSum1+countV;
                walkingTime1=[walkingTime1,walkingCrawlTime(1)];
                crawlingTime1=[crawlingTime1,walkingCrawlTime(2)];
            end           
        end
    end  
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                dataList=[];
                for iTrace=1:size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment,2);
                    dataList=[dataList;[bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(1).Centroid,bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(1).MajorAxisLength]];
                end
                [traceResult,countV,walkingCrawlTime]=traceExploring(dataList,timeInte,velocityBin,walkingLengthThreshold);
                bioTree{iframe}.node{iNode}.Out{iOut}.traceResult=traceResult;
                bioTree{iframe}.node{iNode}.Out{iOut}.velocityDistibution=[velocityBin',countV];
                meanDispTime1(:,1)=meanDispTime1(:,1)+traceResult(:,2);
                meanDispTime1(:,2)=meanDispTime1(:,2)+traceResult(:,3);
                countSum1=countSum1+countV;
                walkingTime1=[walkingTime1,walkingCrawlTime(1)];
                crawlingTime1=[crawlingTime1,walkingCrawlTime(2)];                
            end
        end
    end
end
traceResult.Type0.MSD=[timeInte',meanDispTime0(:,1)./meanDispTime0(:,2)];
traceResult.Type1.MSD=[timeInte',meanDispTime1(:,1)./meanDispTime1(:,2)];
traceResult.Type0.velocity=[velocityBin',countSum0];
traceResult.Type1.velocity=[velocityBin',countSum1];
traceResult.Type0.walkCrawl=[sum(walkingTime0),sum(crawlingTime0)];
traceResult.Type1.walkCrwal=[sum(walkingTime1),sum(crawlingTime1)];
end
function [traceResult,countV,walkingCrawlTime]=traceExploring(dataList,timeInte,velocityBin,walkingLengthThreshold)
velocityCut=10;
meanDisp=zeros(size(timeInte,2),2);
walkingCrawlTime=zeros(1,2);
for iList=1:size(dataList,1)
 if dataList(iList,3)<=walkingLengthThreshold
     walkingCrawlTime(1)=walkingCrawlTime(1)+1;
 else
     walkingCrawlTime(2)=walkingCrawlTime(2)+1;
 end
end
for iInt=1:size(timeInte,2)
    if size(dataList,1)-timeInte(iInt)<=0
        meanDisp(iInt,:)=0;
    end
    if size(dataList,1)-timeInte(iInt)>0
        tempX1=dataList(1:end-timeInte(iInt),1);
        tempY1=dataList(1:end-timeInte(iInt),2);
        tempX2=dataList(1+timeInte(iInt):end,1);
        tempY2=dataList(1+timeInte(iInt):end,2);
        meanDisp(iInt,1)=mean((tempX1-tempX2).^2+(tempY1-tempY2).^2);
        meanDisp(iInt,2)=1;
    end   
end
if size(dataList,1)-velocityCut>1
    tempX1=dataList(1:end-velocityCut,1);
    tempY1=dataList(1:end-velocityCut,2);
    tempX2=dataList(1+velocityCut:end,1);
    tempY2=dataList(1+velocityCut:end,2);
    countV=(histc(((tempX1-tempX2).^2+(tempY1-tempY2).^2),velocityBin));
else
    countV=zeros(size(velocityBin,2),1);
end
traceResult=[timeInte',meanDisp];
end 