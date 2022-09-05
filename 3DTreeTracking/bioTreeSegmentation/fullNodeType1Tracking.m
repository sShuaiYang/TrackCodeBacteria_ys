function [newTrace,usefulNode]=fullNodeType1Tracking(bioTree,nodeInfo,strength,finalOneToOneMatchThreShold)%further track the N input and one Output node, return N new complete bacterial trace
usefulNode=0;
imageSize=bioTree{1}.imageSize;
pixelIdxListIn=getInputMask(bioTree,nodeInfo);
for iIn=1:size(pixelIdxListIn,2)
    newTrace{iIn}.traceInfo.pixelIdxList=[];
    newTrace{iIn}.is2Node=[];
    newTrace{iIn}.isBreakNode=[];
end
for iTrace=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.pixelIdxList,2)
    [pixelIdxListNew,canDivideorNot]=basicDivideSolution(pixelIdxListIn,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.pixelIdxList{iTrace},imageSize,strength,finalOneToOneMatchThreShold);
    if canDivideorNot==1
        leafNum=size(bioTree{iTrace+nodeInfo(1)-2}.leavies,2);
        [newTrace,pixelIdxListIn,isBreak,error]=getNewTrace(newTrace,pixelIdxListNew,iTrace+nodeInfo(1)-2,leafNum);
        if error==1
            usefulNode=1;
            newTrace=[];
            return
        end
        if isBreak==true
            for iIn=1:size(newTrace,2)
                if  isempty(newTrace{iIn}.is2Node)
                    if iTrace<size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.pixelIdxList,2)
                        newTrace{iIn}.traceInfo.pixelIdxList=[newTrace{iIn}.traceInfo.pixelIdxList,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.pixelIdxList(iTrace+1:end)];
                    end
                end
            end
            break;
        end
    else
        if canDivideorNot==0
            if iTrace==1
                usefulNode=1;
                return
            else
                for iIn=1:size(newTrace,2)
                    if isempty(newTrace{iIn}.is2Node)
                        newTrace{iIn}.isBreakNode=1;
                        newTrace{iIn}.breakNodeInfo=[nodeInfo(1)+iTrace-1,0];
                    end
                end
                break;
            end
        end
    end
end
testis2Node=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.is2Node;
if testis2Node==true
    testnodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.nodeInfo;
    inNum=size(bioTree{testnodeInfo(1)}.node{testnodeInfo(2)}.In,2);
end
realInNode=0;
realInLeaf=0;
for iIn=1:size(newTrace,2)
    if ~isempty(newTrace{iIn}.isBreakNode)
        if isempty(newTrace{iIn}.is2Node)
            newTrace{iIn}.is2Node=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.is2Node;
            newTrace{iIn}.afterBreakTrace.pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.pixelIdxList(iTrace:end);
            if newTrace{iIn}.is2Node==1
                newTrace{iIn}.nodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.nodeInfo;
            else
                if newTrace{iIn}.is2Node==0
                    newTrace{iIn}.leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.leafInfo;
                end
            end
        else
%             if isempty(newTrace{iIn}.traceInfo.pixelIdxList)
%                 pixelIdxListIn=getInputMask(bioTree,nodeInfo);
%                 newTrace{iIn}.traceInfo.pixelIdxList{1}=pixelIdxListIn{1,iIn};
%             end
        end
    else
%         if newTrace{iIn}.is2Node==0
%             if isempty(newTrace{iIn}.traceInfo.pixelIdxList)
%                 pixelIdxListIn=getInputMask(bioTree,nodeInfo);
%                 newTrace{iIn}.traceInfo.pixelIdxList{1}=pixelIdxListIn{1,iIn};
%             end
%         end
        if isempty(newTrace{iIn}.is2Node)
            newTrace{iIn}.is2Node=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.is2Node;
            if newTrace{iIn}.is2Node==true
                newTrace{iIn}.nodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.nodeInfo;
                realInNode=realInNode+1;
                if realInNode~=1
                    newTrace{iIn}.nodeInfo(1,3)=inNum+realInNode-1;
%                     newNodeInfo=newTrace{iIn}.nodeInfo;
                end
            end
            if newTrace{iIn}.is2Node==false
                newTrace{iIn}.leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.leafInfo;            
                if realInLeaf~=0
                    newTrace{iIn}.leafInfo(1,2)=size(bioTree{newTrace{iIn}.leafInfo(1)}.leavies,2)+realInLeaf;
                end
                realInLeaf=realInLeaf+1;
            end
        end
    end
end
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
function [newTrace,pixelIdxListIn,isBreak,error]=getNewTrace(newTrace,pixelIdxListNew,iframe,leafNum)
emptyCount=0;
nleaf=0;
for iIn=1:size(newTrace,2)
    if ~isempty(pixelIdxListNew{iIn});
        if ~isempty(newTrace{iIn}.is2Node)
            error=1;
            isBreak=[];
            pixelIdxListIn=[];
            return
        end
        newTrace{iIn}.traceInfo.pixelIdxList=[newTrace{iIn}.traceInfo.pixelIdxList,pixelIdxListNew(iIn)];
    else
        if isempty(newTrace{iIn}.is2Node)
        nleaf=nleaf+1;
        newTrace{iIn}.is2Node=false;
        newTrace{iIn}.leafInfo=[iframe,leafNum+nleaf];
        end
        emptyCount= emptyCount+1;
    end
end
if size(newTrace,2)-emptyCount==1
    isBreak=true;
else
    isBreak=false;
end
error=0;
pixelIdxListIn=pixelIdxListNew;
end
