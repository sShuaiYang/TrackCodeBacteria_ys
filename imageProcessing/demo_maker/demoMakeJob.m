function demoImage=demoMakeJob(traceInfo,segreationImage,maskImage)
startFrame=1;
endFrame=1665;
divisionColor=[0,0,1];
attachingColor=[1,0,0];
detchingColor=[0,1,0];

demoImage=zeros(size(segreationImage,1),size(segreationImage,2),3,endFrame-startFrame+1,'uint8');

parfor iframe=startFrame:endFrame
    maskImageinFrame=maskImage(:,:,iframe);
    traceInframe=traceInfo([traceInfo.startFrame]<=iframe & [traceInfo.endFrame]>=iframe);
    infoInframe=getInfoinframe(traceInframe,iframe);
    xPos=infoInframe(1,:);
    yPos=infoInframe(2,:);
    
    findDivisonindxe=[infoInframe(6,:)]==-2;
    xDivision=xPos(findDivisonindxe);
    yDivision=yPos(findDivisonindxe);
    maskDivision=bwselect(maskImageinFrame,xDivision,yDivision,4);
    demoImageTemp=imoverlay(segreationImage(:,:,iframe),maskDivision,divisionColor);
    
    findAttachingindxe=[infoInframe(6,:)]==-1;
    xAttaching=xPos(findAttachingindxe);
    yAttaching=yPos(findAttachingindxe);
    maskAttaching=bwselect(maskImageinFrame,xAttaching,yAttaching,4);
    demoImageTemp=imoverlay(demoImageTemp,maskAttaching,attachingColor);
    
    findDetchingindxe=[infoInframe(5,:)]==-1;
    endTraceframe=infoInframe(4,:);
    findDetchingindxe=endTraceframe(findDetchingindxe)==iframe;
    xDetching=xPos(findDetchingindxe);
    yDetching=yPos(findDetchingindxe);
    maskDetching=bwselect(maskImageinFrame,xDetching,yDetching,4);
    demoImage(:,:,:,iframe)=imoverlay(demoImageTemp,maskDetching,detchingColor);
    
    
end
end

function infoInframe=getInfoinframe(traceInframe,iframe)
beadNumberinframe=length(traceInframe);
infoInframe=zeros(6,beadNumberinframe);
startTraceframe=[traceInframe.startFrame];
endTraceframe=[traceInframe.endFrame];
possibleChildren=[traceInframe.possibleChildren];
possibleParent=[traceInframe.possibleParent];
for iBead=1:beadNumberinframe
    startBeadframe=traceInframe(iBead).startFrame;
    xyPos=traceInframe(iBead).frameDetail(iframe-startBeadframe+1).Centroid;
    xPos=xyPos(1);
    yPos=xyPos(2);
    infoInframe(1,iBead)=xPos;
    infoInframe(2,iBead)=yPos;
end
infoInframe(3,:)=startTraceframe;
infoInframe(4,:)=endTraceframe;
infoInframe(5,:)=possibleChildren;
infoInframe(6,:)=possibleParent;
end
