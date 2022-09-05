function distanceResult=spatialDistance(allData,paraName,type)
if strcmp(type,'diff')
    % MSDslope,retentionTime,averageLength,aveGrowthRate,aveVelocity
    distanceResult=[];
    for iNum=1:numel(allData)
        dtInfo=allData(iNum);
        dataNeed=false(numel(dtInfo.pieceDetail));
        for i=1:numel(dtInfo.pieceDetail)
            if dtInfo.pieceDetail(i).retentionTime>=1000
                dataNeed(i)=1;
            end
        end
        focusInfo=dtInfo.pieceLinkInfo(dataNeed,:);
        orderNum=1:size(dtInfo.linkMatrix,1);
        pieceDetail=dtInfo.pieceDetail(dataNeed);
        for i=1:size(dtInfo.linkMatrix,1)
            for j=i+1:size(dtInfo.linkMatrix,1)
                if dtInfo.linkMatrix(i,j)==1
                    dtInfo.linkMatrix(j,i)=1;
                end
            end
        end
        dtInfo.linkMatrix=sparse(dtInfo.linkMatrix);
        for i=1:size(focusInfo,1)
            for j=i+1:size(focusInfo,1)
                headOrder=dtInfo.allList(:,1)==focusInfo(i,1) & dtInfo.allList(:,2)==focusInfo(i,2) & dtInfo.allList(:,3)==focusInfo(i,3);
                tailOrder=dtInfo.allList(:,1)==focusInfo(j,4) & dtInfo.allList(:,2)==focusInfo(j,5) & dtInfo.allList(:,3)==focusInfo(j,6);
                [dist,~,~]=graphshortestpath(dtInfo.linkMatrix,orderNum(headOrder),orderNum(tailOrder));
                if dist==inf
                    p=1;
                end
                eval(strcat('diffResult=abs(pieceDetail(i).',paraName,'-pieceDetail(j).',paraName,');'))
                distanceResult=[distanceResult;dist,diffResult];
            end
        end
    end
    plot(distanceResult(:,1),distanceResult(:,2),'lineStyle','none','marker','o');
    dis=distanceResult(:,1);
    disResult=distanceResult(:,2);
    for i=1:max(dis)
        result(i)=mean(disResult(dis==i));
    end
    hold on;plot(result,'r')
end
if strcmp(type,'normal')
    distanceResult=[];
%         for iNum=1:numel(allData)
    for iNum=2
        dtInfo=allData(iNum);
        dataNeed=false(numel(dtInfo.pieceDetail));
        for i=1:numel(dtInfo.pieceDetail)
            if dtInfo.pieceDetail(i).retentionTime>=1000
                dataNeed(i)=1;
            end
        end
        pieceDetail=dtInfo.pieceDetail(dataNeed);
        edgeGeneral=dtInfo.edgeGeneral(dataNeed);
        for i=1:numel(pieceDetail)
            eval(strcat('aimPara=pieceDetail(i).',paraName,';'));
            distanceResult=[distanceResult;edgeGeneral(i),aimPara];
        end
    end
    plot(distanceResult(:,1),distanceResult(:,2),'lineStyle','none','marker','o');
    dis=distanceResult(:,1);
    disResult=distanceResult(:,2);
    for i=1:max(dis)
        result(i)=mean(disResult(dis==i));
    end
    hold on;plot(result,'r')
end
end
        
