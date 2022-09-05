function bioTree1=makeBioTree(treeBranch,startFrame,endFrame,frameShift,xSize,ySize)
fprintf('\n');
bioTree1=twoFrameConnect(treeBranch(startFrame),startFrame+frameShift);
for iConnect=startFrame+1:endFrame-1
    dispFrame(iConnect+frameShift);
    bioTree2=twoFrameConnect(treeBranch(iConnect),iConnect+frameShift);
    if ~isempty(bioTree1{iConnect+frameShift}.leavies)
        bioTree1=bioTreeConnecter(bioTree1,bioTree2,iConnect+frameShift,xSize,ySize);
    else
        bioTree2(1:iConnect+frameShift)=bioTree1;
        bioTree1=bioTree2;
    end
end

end

function dispFrame(iConnect) 
lengthB=length(num2str(iConnect-1));
back='';
for i=1:lengthB
    back=strcat(back,'\b');
end
fprintf(strcat(back,num2str(iConnect)));
end