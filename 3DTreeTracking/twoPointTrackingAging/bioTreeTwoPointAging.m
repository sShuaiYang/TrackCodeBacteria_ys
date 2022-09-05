function bioTree=bioTreeTwoPointAging(bioTree)
bioTree=bioTreeTwoPointTracking(bioTree);
bioTree=bioTreePointMatch(bioTree);
end
function bioTree=bioTreeTwoPointTracking(bioTree)
bioTree=getTwoPointPosition(bioTree);
end
function bioTree=getTwoPointPosition(bioTree)
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if ~isempty(bioTree{iframe}.root{iroot}.traceInfo.measurment)
                for iMeasurment=1:size(bioTree{iframe}.root{iroot}.traceInfo.measurment{1},1)
                    centroid=bioTree{iframe}.root{iroot}.traceInfo.measurment{1}(iMeasurment).Centroid;
                    orientation=bioTree{iframe}.root{iroot}.traceInfo.measurment{1}(iMeasurment).Orientation;
                    majorAxisLength=bioTree{iframe}.root{iroot}.traceInfo.measurment{1}(iMeasurment).MajorAxisLength;
                    eccentricity=bioTree{iframe}.root{iroot}.traceInfo.measurment{1}(iMeasurment).Eccentricity;
                    [p1Position,p2Position]=twoPointPosition(centroid,orientation,majorAxisLength,eccentricity);
                    bioTree{iframe}.root{iroot}.traceInfo.measurment{1}(iMeasurment).p1Position=p1Position;
                    bioTree{iframe}.root{iroot}.traceInfo.measurment{1}(iMeasurment).p2Position=p2Position;
                end
                for iTrace=2:size(bioTree{iframe}.root{iroot}.traceInfo.measurment,2)
                    for iMeasurment=1:size(bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace},1)
                        centroid_next=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(iMeasurment).Centroid;
                        orientation_next=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(iMeasurment).Orientation;
                        majorAxisLength_next=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(iMeasurment).MajorAxisLength;
                        eccentricity_next=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(iMeasurment).Eccentricity;
                        [p1Position_next,p2Position_next]=twoPointPosition(centroid_next,orientation_next,majorAxisLength_next,eccentricity_next);
                        p1Position_pre= bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace-1}(1).p1Position;
                        p2Position_pre= bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace-1}(1).p2Position;
                        [p1Position,p2Position]=twoPointConnect(p1Position_pre,p1Position_next,p2Position_pre,p2Position_next);
                        bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(iMeasurment).p1Position=p1Position;
                        bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(iMeasurment).p2Position=p2Position;
                    end
                end
            end
        end
    end    
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                if ~isempty(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment)
                    for iMeasurment=1:size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1},1)
                        centroid= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}(iMeasurment).Centroid;
                        orientation= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}(iMeasurment).Orientation;
                        majorAxisLength= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}(iMeasurment).MajorAxisLength;
                        eccentricity= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}(iMeasurment).Eccentricity;
                        [p1Position,p2Position]=twoPointPosition(centroid,orientation,majorAxisLength,eccentricity);
                        bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}(iMeasurment).p1Position=p1Position;
                        bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}(iMeasurment).p2Position=p2Position;
                    end
                    for iTrace=2:size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment,2)
                        for iMeasurment=1:size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace},1)
                            centroid_next= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(iMeasurment).Centroid;
                            orientation_next= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(iMeasurment).Orientation;
                            majorAxisLength_next= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(iMeasurment).MajorAxisLength;
                            eccentricity_next= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(iMeasurment).Eccentricity;
                            [p1Position_next,p2Position_next]=twoPointPosition(centroid_next,orientation_next,majorAxisLength_next,eccentricity_next);
                            p1Position_pre=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace-1}(1).p1Position;
                            p2Position_pre=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace-1}(1).p2Position;
                            [p1Position,p2Position]=twoPointConnect(p1Position_pre,p1Position_next,p2Position_pre,p2Position_next);
                            bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(iMeasurment).p1Position=p1Position;
                            bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(iMeasurment).p2Position=p2Position;
                        end
                    end
                end
            end
        end
    end
end
end
function [p1Position,p2Position]=twoPointPosition(centroid,orientation,majorAxisLength,eccentricity)
O_sign=sign(orientation);
if orientation==0
    O_sign=1;
end
O_cos=cos(abs(pi*orientation/180));
O_sin=sin(abs(pi*orientation/180));
Axis_length=0.5*majorAxisLength;
x=centroid(1);
y=centroid(2);
p1_x=x-O_sign*Axis_length*O_cos*eccentricity;
p1_y=y+Axis_length*O_sin*eccentricity;
p2_x=x+O_sign*Axis_length*O_cos*eccentricity;
p2_y=y-Axis_length*O_sin*eccentricity;
p1Position=[p1_x,p1_y];
p2Position=[p2_x,p2_y];
end
function  [p1Position,p2Position]=twoPointConnect(p1Position_pre,p1Position_next,p2Position_pre,p2Position_next)
vectorP1P2_pre=p1Position_pre-p2Position_pre;
vectorP1P2_next=p1Position_next-p2Position_next;
if sum(vectorP1P2_pre.*vectorP1P2_next)>0
    p1Position=p1Position_next;
    p2Position=p2Position_next;
    return;
end
if sum(vectorP1P2_pre.*vectorP1P2_next)<0
    p1Position=p2Position_next;
    p2Position=p1Position_next;
    return;
end
if sum(vectorP1P2_pre.*vectorP1P2_next)==0
    D1=sum((p1Position_pre-p1Position_next).^2);
    D2=sum((p1Position_pre-p2Position_next).^2);
    if D1<=D2
        p1Position=p1Position_next;
        p2Position=p2Position_next;
        return;
    else
        p1Position=p2Position_next;
        p2Position=p1Position_next;
        return;
    end
end
end
function bioTree=bioTreePointMatch(bioTree)
for iframe=1:size(bioTree,2)
    for iRoot=1:size(bioTree{iframe}.root,2)
        bioTree{iframe}.root{iRoot}.isOld1=1;
        bioTree{iframe}.root{iRoot}.isOld2=1;
    end
end
for iframe=1:size(bioTree,2)
    for iNode=1:size(bioTree{iframe}.node,2)
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)            
            bioTree{iframe}.node{iNode}.Out{iOut}.isOld1=-inf;
            bioTree{iframe}.node{iNode}.Out{iOut}.isOld2=-inf;
        end
    end
end
for iframe=1:size(bioTree,2)
    for iNode=1:size(bioTree{iframe}.node,2)
        if size(bioTree{iframe}.node{iNode}.In,2)==1 && size(bioTree{iframe}.node{iNode}.Out,2)~=1
            preIsNode=bioTree{iframe}.node{iNode}.In{1}.isNode;
            if preIsNode==0
                outInfo=[];
                preRoot=bioTree{iframe}.node{iNode}.In{1}.rootInfo;
                inInfo=[bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.measurment{end}.p1Position,bioTree{preRoot(1)}.root{preRoot(2)}.isOld1;...
                    bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.measurment{end}.p2Position,bioTree{preRoot(1)}.root{preRoot(2)}.isOld2];
                for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                    outInfo(iOut*2-1:2*iOut,:)=[bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}.p1Position;bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}.p2Position];
                end
                distMatrix=pdist2(inInfo(:,1:2),outInfo);
                orderNum=1:size(outInfo,1);
                matchNum=orderNum(distMatrix(1,:)==min(distMatrix(1,:)));
                matchNum=matchNum(1);
                outNum=fix((matchNum+1)/2);
                leaveNum=matchNum-(outNum-1)*2;
                if leaveNum==1
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=bioTree{preRoot(1)}.root{preRoot(2)}.isOld1+1;
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=0;
                end
                if leaveNum==2
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=0;
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=bioTree{preRoot(1)}.root{preRoot(2)}.isOld1+1;
                end
                distMatrix(:,matchNum)=50;
                matchNum=orderNum(distMatrix(2,:)==min(distMatrix(2,:)));
                matchNum=matchNum(1);
                outNum=fix((matchNum+1)/2);
                leaveNum=matchNum-(outNum-1)*2;
                if leaveNum==1
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=bioTree{preRoot(1)}.root{preRoot(2)}.isOld2+1;
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=0;
                end
                if leaveNum==2
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=0;
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=bioTree{preRoot(1)}.root{preRoot(2)}.isOld2+1;
                end
            end
            if preIsNode==1
                outInfo=[];
                preNode=bioTree{iframe}.node{iNode}.In{1}.nodeInfo;
                inInfo=[bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.measurment{end}.p1Position;bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.measurment{end}.p2Position];
                for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                    outInfo(iOut*2-1:2*iOut,:)=[bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}.p1Position;bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}.p2Position];
                end
                distMatrix=pdist2(inInfo(:,1:2),outInfo);
                orderNum=1:size(outInfo,1);
                if ~isempty(bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.isOld1)
                    matchNum=orderNum(distMatrix(1,:)==min(distMatrix(1,:)));
                    matchNum=matchNum(1);
                    outNum=fix((matchNum+1)/2);
                    leaveNum=matchNum-(outNum-1)*2;
                    if leaveNum==1
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.isOld1+1;
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=0;
                    end
                    if leaveNum==2
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=0;
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.isOld1+1;
                    end
                    distMatrix(:,matchNum)=50;
                    matchNum=orderNum(distMatrix(2,:)==min(distMatrix(2,:)));
                    matchNum=matchNum(1);
                    outNum=fix((matchNum+1)/2);
                    leaveNum=matchNum-(outNum-1)*2;
                    if leaveNum==1
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.isOld2+1;
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=0;
                    end
                    if leaveNum==2
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=0;
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.isOld2+1;
                    end
                end
            end
        end
    end
end
end