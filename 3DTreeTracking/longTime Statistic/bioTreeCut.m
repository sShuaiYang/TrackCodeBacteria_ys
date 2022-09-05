function bioTree=bioTreeCut(bioTree,myFrame)
for iFrame=1:myFrame
    for iRoot=1:size(bioTree{iFrame}.root,2)
        endLeafNum=size(bioTree{myFrame}.leavies,2);
        if bioTree{iFrame}.root{iRoot}.is2Node==1
            nodeInfo=bioTree{iFrame}.root{iRoot}.nodeInfo;
            if nodeInfo(1)>myFrame
                bioTree{iFrame}.root{iRoot}.is2Node=0;
                bioTree{iFrame}.root{iRoot}.nodeInfo=[];
                leafInfo=[myFrame,size(bioTree{myFrame}.leavies,2)+1];
                bioTree{iFrame}.root{iRoot}.leafInfo=leafInfo;
                deadLineNum=myFrame-iFrame+1;
                bioTree{iFrame}.root{iRoot}.traceInfo.pixelIdxList(deadLineNum+1:end)=[];
%                 bioTree{iFrame}.root{iRoot}.traceInfo.measurment(deadLineNum+1:end)=[];
                bioTree{myFrame}.leavies{endLeafNum+1}.is2Node=0;
                bioTree{myFrame}.leavies{endLeafNum+1}.nodeInfo=[];
                bioTree{myFrame}.leavies{endLeafNum+1}.rootInfo=[iFrame,iRoot];
                bioTree{myFrame}.leavies{endLeafNum+1}.leaviesPixelDetail=bioTree{iFrame}.root{iRoot}.traceInfo.pixelIdxList{end};
%                 bioTree{myFrame}.leavies{endLeafNum+1}.leafMeasurment=bioTree{iFrame}.root{iRoot}.traceInfo.measurment{end};
            end
        end
        if bioTree{iFrame}.root{iRoot}.is2Node==0
            leafInfo=bioTree{iFrame}.root{iRoot}.leafInfo;
            if leafInfo(1)>myFrame
                bioTree{iFrame}.root{iRoot}.is2Node=0;
                bioTree{iFrame}.root{iRoot}.nodeInfo=[];
                leafInfo=[myFrame,size(bioTree{myFrame}.leavies,2)+1];
                bioTree{iFrame}.root{iRoot}.leafInfo=leafInfo;
                deadLineNum=myFrame-iFrame+1;
                bioTree{iFrame}.root{iRoot}.traceInfo.pixelIdxList(deadLineNum+1:end)=[];
%                 bioTree{iFrame}.root{iRoot}.traceInfo.measurment(deadLineNum+1:end)=[];
                bioTree{myFrame}.leavies{endLeafNum+1}.is2Node=0;
                bioTree{myFrame}.leavies{endLeafNum+1}.nodeInfo=[];
                bioTree{myFrame}.leavies{endLeafNum+1}.rootInfo=[iFrame,iRoot];
                bioTree{myFrame}.leavies{endLeafNum+1}.leaviesPixelDetail=bioTree{iFrame}.root{iRoot}.traceInfo.pixelIdxList{end};
%                 bioTree{myFrame}.leavies{endLeafNum+1}.leafMeasurment=bioTree{iFrame}.root{iRoot}.traceInfo.measurment{end};
            end
        end
    end
    for iNode=1:size(bioTree{iFrame}.node,2)
        for iOut=1:size(bioTree{iFrame}.node{iNode}.Out,2)
            endLeafNum=size(bioTree{myFrame}.leavies,2);
            if bioTree{iFrame}.node{iNode}.Out{iOut}.is2Node==1
                nodeInfo=bioTree{iFrame}.node{iNode}.Out{iOut}.nodeInfo;
                %                 if nodeInfo(1)>myFrame
                if nodeInfo(1)>myFrame&&~isempty(bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList)%ys 20191230
                    bioTree{iFrame}.node{iNode}.Out{iOut}.is2Node=0;
                    bioTree{iFrame}.node{iNode}.Out{iOut}.nodeInfo=[];
                    leafInfo=[myFrame,size(bioTree{myFrame}.leavies,2)+1];
                    bioTree{iFrame}.node{iNode}.Out{iOut}.leafInfo=leafInfo;
                    deadLineNum=myFrame-iFrame+1;
                    bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList(deadLineNum+1:end)=[];
                    %                     bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo.measurment(deadLineNum+1:end)=[];
                    bioTree{myFrame}.leavies{endLeafNum+1}.is2Node=1;
                    bioTree{myFrame}.leavies{endLeafNum+1}.rootInfo=[];
                    bioTree{myFrame}.leavies{endLeafNum+1}.nodeInfo=[iFrame,iNode,iOut];
                    
                    bioTree{myFrame}.leavies{endLeafNum+1}.leaviesPixelDetail=bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList{end};
                    %                     bioTree{myFrame}.leavies{endLeafNum+1}.leafMeasurment=bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo.measurment{end};
                end
            end
            if bioTree{iFrame}.node{iNode}.Out{iOut}.is2Node==0
                leafInfo=bioTree{iFrame}.node{iNode}.Out{iOut}.leafInfo;
                if leafInfo(1)>myFrame
                    bioTree{iFrame}.node{iNode}.Out{iOut}.is2Node=0;
                    bioTree{iFrame}.node{iNode}.Out{iOut}.nodeInfo=[];
                    leafInfo=[myFrame,size(bioTree{myFrame}.leavies,2)+1];
                    bioTree{iFrame}.node{iNode}.Out{iOut}.leafInfo=leafInfo;
                    deadLineNum=myFrame-iFrame+1;
                    bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList(deadLineNum+1:end)=[];
%                     bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo.measurment(deadLineNum+1:end)=[];
                    bioTree{myFrame}.leavies{endLeafNum+1}.is2Node=1;
                    bioTree{myFrame}.leavies{endLeafNum+1}.rootInfo=[];
                    bioTree{myFrame}.leavies{endLeafNum+1}.nodeInfo=[iFrame,iNode,iOut];
                    bioTree{myFrame}.leavies{endLeafNum+1}.leaviesPixelDetail=bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList{end};
%                     bioTree{myFrame}.leavies{endLeafNum+1}.leafMeasurment=bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo.measurment{end};
                end
            end
        end
    end
end
bioTree(myFrame+1:end)=[];
end
            