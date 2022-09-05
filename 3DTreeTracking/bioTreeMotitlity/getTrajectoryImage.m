function trajectoryImage=getTrajectoryImage(bioTree,startFrame,stepFrame,endFrame)
xSize=bioTree{1}.imageSize(1);
ySize=bioTree{1}.imageSize(2);
timeStack=startFrame:stepFrame:endFrame;
timeSize=size(timeStack,2);
imageP1=zeros(xSize,ySize,timeSize,'uint8');
imageP2=zeros(xSize,ySize,timeSize,'uint8');
maskImage=false(xSize,ySize,timeSize);

% imageCentroid=zeros(xSize,ySize,timeSize,'int16');
for iframe=1:endFrame
    if iframe >endFrame
        break;
    end
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            for iTrace=1:size(bioTree{iframe}.root{iroot}.traceInfo.measurment,2)
                timePoint=getTimePoint(iframe,iTrace,startFrame,stepFrame,endFrame);
                if ~isempty(timePoint)
                    if timeStack(timePoint)==iframe+iTrace-1
                        tempImage=false(xSize,ySize);
                        pixelIdxList=bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList{iTrace};
                        tempImage(pixelIdxList)=true;
                        maskImage(:,:,timePoint)=maskImage(:,:,timePoint)|tempImage;
                    end
                    %                 centroid=round(bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(1).Centroid);
                    p1Position=round(bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(1).p1Position);
                    p2Position=round(bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(1).p2Position);
                    %                     imageCentroid(centroid(1),centroid(2),timePoint)=imageCentroid(centroid(1),centroid(2))+1;
                    if p1Position(1)>=0 &&p1Position(1)<ySize && p1Position(2)>=0 && p1Position(2)<xSize
                        imageP1(p1Position(2),p1Position(1),timePoint)=imageP1(p1Position(2),p1Position(1),timePoint)+1;
                    end
                    if p2Position(1)>=0 &&p2Position(1)<ySize && p2Position(2)>=0 && p2Position(2)<xSize
                        imageP2(p2Position(2),p2Position(1),timePoint)=imageP2(p2Position(2),p2Position(1),timePoint)+1;
                    end
                end
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                for iTrace=1:size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment,2)
                    timePoint=getTimePoint(iframe,iTrace,startFrame,stepFrame,endFrame);
                    if ~isempty(timePoint)
                        if timeStack(timePoint)==iframe+iTrace-1
                            tempImage=false(xSize,ySize);
                            pixelIdxList=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList{iTrace};
                            tempImage(pixelIdxList)=true;
                            maskImage(:,:,timePoint)=maskImage(:,:,timePoint)|tempImage;
                        end
                        %                 centroid=round(bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(1).Centroid);
                        p1Position=round(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(1).p1Position);
                        p2Position=round(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(1).p2Position);
                        %                     imageCentroid(centroid(1),centroid(2),timePoint)=imageCentroid(centroid(1),centroid(2))+1;
                        if p1Position(1)>=0 &&p1Position(1)<ySize && p1Position(2)>=0 && p1Position(2)<xSize
                            imageP1(p1Position(2),p1Position(1),timePoint)=imageP1(p1Position(2),p1Position(1),timePoint)+1;
                        end
                        if p2Position(1)>=0 &&p2Position(1)<ySize && p2Position(2)>=0 && p2Position(2)<xSize
                            imageP2(p2Position(2),p2Position(1),timePoint)=imageP2(p2Position(2),p2Position(1),timePoint)+1;
                        end
                    end
                end
            end
        end
    end
end
trajectoryImage=ovelayImage( maskImage,imageP1,imageP2);
end
function timePoint=getTimePoint(iframe,iTrace,startFrame,stepFrame,endFrame)
realTime=iframe+iTrace-1;
timeFrame=startFrame:stepFrame:endFrame;
timePoint=timeFrame(timeFrame>realTime-stepFrame & timeFrame <= realTime);
timePoint=find(timeFrame==timePoint);
end
function trajectoryImage=ovelayImage( maskImage,imageP1,imageP2)
for iframe=1:size(maskImage,3)
    maskImage(:,:,iframe)=imfill(maskImage(:,:,iframe),'holes');
    maskImage(:,:,iframe)=bwmorph(maskImage(:,:,iframe),'dilate');
    maskImage(:,:,iframe)=bwmorph(maskImage(:,:,iframe),'remove');
end
trajectoryImage=zeros(size(maskImage,1),size(maskImage,2),3,size(maskImage,3),'single');
trajectoryImage(:,:,1,1)=(imageP1(:,:,1));
trajectoryImage(:,:,2,1)=(imageP2(:,:,1));
trajectoryImage(:,:,3,1)=maskImage(:,:,1);
for iframe=2:size(maskImage,3)
    trajectoryImage(:,:,1,iframe)=(single(imageP1(:,:,iframe)))+ trajectoryImage(:,:,1,iframe-1);
    trajectoryImage(:,:,2,iframe)=(single(imageP2(:,:,iframe)))+ trajectoryImage(:,:,2,iframe-1);
    trajectoryImage(:,:,3,iframe)=single(maskImage(:,:,iframe));
end
for iframe=1:size(maskImage,3)
    trajectoryImage(:,:,1,iframe)=log(trajectoryImage(:,:,1,iframe)+1);
    trajectoryImage(:,:,2,iframe)=log(trajectoryImage(:,:,2,iframe)+1);
    sumP1=max(max(trajectoryImage(:,:,1,iframe)));
    sumP2=max(max(trajectoryImage(:,:,2,iframe)));
    trajectoryImage(:,:,1,iframe)=  trajectoryImage(:,:,1,iframe)./sumP1;
    trajectoryImage(:,:,2,iframe)=  trajectoryImage(:,:,2,iframe)./sumP2;
end
end