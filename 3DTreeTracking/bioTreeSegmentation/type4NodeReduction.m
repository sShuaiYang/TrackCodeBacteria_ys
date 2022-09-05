function bioTree=type4NodeReduction(bioTree)
nodeCount=0;
for iframe=1:size(bioTree,2)
% disp(strcat(num2str(iframe),'_'))
% fprintf('\b\n')
dispFrame(iframe)
if ~isempty(bioTree{iframe}.node)
    nodeList=type4NodeinFrame(bioTree{iframe},iframe);
    if ~isempty(nodeList)
        for iList=1:size(nodeList,1)
            %nodeInfoIn=getInputInfo(bioTree,nodeList(iList,:));
            nodeInfo=nodeList(iList,:);
            pixelIdxListIn=getInputMask(bioTree,nodeInfo);
            iInProps=0;
            for iIn=1:size(pixelIdxListIn,2)
                [~,BWImage]=idx2Xy(pixelIdxListIn{iIn},bioTree{1}.imageSize);
                CC=bwconncomp(BWImage);
                regionNum=CC.NumObjects;
                if regionNum>=2
                    iInProps=1;
                end
            end
            if iInProps==1
                break
            end
            imageSize=bioTree{1}.imageSize;
            traceInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.pixelIdxList;
            pictureOut=false(imageSize);
            pictureOut(traceInfo{1})=true;
            cc=bwconncomp(pictureOut);
            if cc.NumObjects==2
                outPixelIdxList=cc.PixelIdxList;
                if numel(traceInfo)>1
                    nodeNum=size(bioTree{nodeInfo(1)+1}.node,2);
                    bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1};
                    bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.traceInfo.pixelIdxList(1)=[];
                    if bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.is2Node==1
                        newNodeInfo=bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.nodeInfo;
                        bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{newNodeInfo(3)}.nodeInfo=[nodeInfo(1)+1,nodeNum+1,1];
                    end
                    if bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.is2Node==0
                        newLeafInfo=bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.leafInfo;
                        bioTree{newLeafInfo(1)}.leavies{newLeafInfo(2)}.nodeInfo=[nodeInfo(1)+1,nodeNum+1,1];
                    end
                    for iOut=1:2
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.is2Node=1;
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.leafInfo=[];
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.nodeInfo=[nodeInfo(1)+1,nodeNum+1,iOut];
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.traceInfo.pixelIdxList=[];
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.traceInfo.pixelIdxList{1}=outPixelIdxList{iOut};
                    end
                    for iIn=1:2
                        bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.In{1,iIn}.isNode=1;
                        bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.In{1,iIn}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                    end
                end
                if numel(traceInfo)==1
                    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.is2Node==1
                        newNodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.nodeInfo;
                        inNum=size(bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In,2);
                        for iIn=1:2
                            if iIn==2
                                bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,iIn+inNum-1}.isNode=1;
                                bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,iIn+inNum-1}.rootInfo=[];
                                bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,iIn+inNum-1}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                            else
                                bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,newNodeInfo(3)}.isNode=1;
                                bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,newNodeInfo(3)}.rootInfo=[];
                                bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,newNodeInfo(3)}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                            end
                        end
                        for iOut=1:2
                            if iOut==2
                                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.is2Node=1;
                                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.leafInfo=[];
                                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.nodeInfo=[newNodeInfo(1),newNodeInfo(2),iOut+inNum-1];
                                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.traceInfo.pixelIdxList{1}=outPixelIdxList{iOut};
                            else
                                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.is2Node=1;
                                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.leafInfo=[];
                                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.nodeInfo=newNodeInfo;
                                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.traceInfo.pixelIdxList{1}=outPixelIdxList{iOut};
                            end
                        end
                    end
                end
            end
            nodeCount=nodeCount+1;
%             disp(nodeCount)
        end
    end
end
end
fprintf('\n')
end
function dispFrame(iConnect) 
lengthB=length(num2str(iConnect-1));
back='';
for i=1:lengthB
    back=strcat(back,'\b');
end
fprintf(strcat(back,num2str(iConnect)));
end
function nodeList=type4NodeinFrame(bioTreeFrame,iframe)
nodeList=[];
for iNode=1:size(bioTreeFrame.node,2)
    if size(bioTreeFrame.node{iNode}.In,2)==1 && size(bioTreeFrame.node{iNode}.Out,2)==1
        nodeList=[nodeList;[iframe,iNode]];
    end
end
end
function [xyMin,BWImageGain]=idx2Xy(pixelIdxList,pictureSize)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=round(pixelIdxList-(yresult-1)*xSize);
xMin=min(xresult);
xMax=max(xresult);
yMin=min(yresult);
yMax=max(yresult);
yresult2=yresult-yMin+1;
xresult2=xresult-xMin+1;
BWImage=false(xMax-xMin+1,yMax-yMin+1);
pixelIdxList2=xresult2+(yresult2-1)*(xMax-xMin+1);
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
xyMin=[xMin,yMin];
BWImageGain(2:end-1,2:end-1)=BWImage;
end
function pixelIdxListIn=getInputMask(bioTree,nodeInfo)
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
        nodeInfo_pre=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
        pixelIdxListIn{iIn}=bioTree{nodeInfo_pre(1)}.node{nodeInfo_pre(2)}.Out{nodeInfo_pre(3)}.traceInfo.pixelIdxList{end};
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
        rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
        if rootInfo(1)==nodeInfo(1)
            pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.rootPixelDetail;
        end
        if  rootInfo(1)<nodeInfo(1)
            pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{end};
        end
    end
end
end