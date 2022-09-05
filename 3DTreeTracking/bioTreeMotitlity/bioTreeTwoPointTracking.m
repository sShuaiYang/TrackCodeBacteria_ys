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