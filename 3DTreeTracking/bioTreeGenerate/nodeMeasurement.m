function bioTree=nodeMeasurement(bioTree,frameShift)
divisionThreshold=0.1;
frameThreshold=150;
bioTree=classNode(bioTree,frameShift);
bioTree=nodeOrder(bioTree,frameShift);
bioTree=markDivison(bioTree,frameShift,divisionThreshold,frameThreshold);
divisionCount(bioTree,frameShift);
end
function bioTree=classNode(bioTree,frameShift)
type1Count=0;
type2Count=0;
type0Count=0;
for iframe=frameShift+1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            bioTree{iframe}.node{iNode}.type=[];
            nodeIn=size(bioTree{iframe}.node{iNode}.In,2);
            nodeOut=size(bioTree{iframe}.node{iNode}.Out,2);
            if nodeIn>nodeOut
                bioTree{iframe}.node{iNode}.type=1;
                type1Count=type1Count+1;
            end
            if nodeIn<nodeOut
                bioTree{iframe}.node{iNode}.type=2;
                type2Count=type2Count+1;
            end
            if nodeIn==nodeOut
                bioTree{iframe}.node{iNode}.type=0;
                type0Count=type0Count+1;
            end
        end
    end
end
disp(strcat('type1Node=',num2str(type1Count)));
disp(strcat('type2Node=',num2str(type2Count)));
disp(strcat('type0Node=',num2str(type0Count)));
end
function bioTree=nodeOrder(bioTree,frameShift)
for iframe=frameShift+1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            orderNum=[];
            parentNode=[];
            for iIn=1:size(bioTree{iframe}.node{iNode}.In,2)
                if bioTree{iframe}.node{iNode}.In{iIn}.isNode==false
                    orderNum=[orderNum,0];
                    parentNode=[parentNode;bioTree{iframe}.node{iNode}.In{iIn}.rootInfo];
                else
                    nodeInfo=bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo;
                    orderNum=[orderNum,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.order+1];
                    parentNode=[parentNode;[bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo(1),bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo(2)]];
                end
            end
            [bioTree{iframe}.node{iNode}.order,idx]=max(orderNum);
            bioTree{iframe}.node{iNode}.parentNode=parentNode(idx,:);
        end
    end
end
end
function bioTree=markDivison(bioTree,frameShift,divisionThreshold,frameThreshold)
for iframe=frameShift+1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            if bioTree{iframe}.node{iNode}.type==2 && (size(bioTree{iframe}.node{iNode}.Out,2)-size(bioTree{iframe}.node{iNode}.In,2))==1
                if bioTree{iframe}.node{iNode}.order==0
                    bioTree{iframe}.node{iNode}.type=3;
                    continue;
                else
                    parentInfo=findParent(bioTree,iframe,iNode);
                    if isempty(parentInfo)
                        fillArea1=bioTree{iframe}.node{iNode}.Out{1}.traceInfo.measurment{1}.FilledArea;
                        fillArea2=bioTree{iframe}.node{iNode}.Out{2}.traceInfo.measurment{1}.FilledArea;
                        if  abs((fillArea1-fillArea2))/(fillArea1+fillArea2)<=divisionThreshold
                            bioTree{iframe}.node{iNode}.type=3;
                            continue;
                        end
                    else
                        fillArea1=bioTree{iframe}.node{iNode}.Out{1}.traceInfo.measurment{1}.FilledArea;
                        fillArea2=bioTree{iframe}.node{iNode}.Out{2}.traceInfo.measurment{1}.FilledArea;
                        if  abs((fillArea1-fillArea2))/(fillArea1+fillArea2)<=divisionThreshold && (iframe-parentInfo(1))>=frameThreshold
                            bioTree{iframe}.node{iNode}.type=3;
                            continue;
                        end
                    end
                end
            end
        end
    end
end
end
function parentInfo=findParent(bioTree,iframe,iNode)
parentInfo=[];
parentTempInfo=bioTree{iframe}.node{iNode}.parentNode;
for iOrder=1:bioTree{iframe}.node{iNode}.order
    if bioTree{parentTempInfo(1)}.node{parentTempInfo(2)}.type==3
        parentInfo=parentTempInfo;
        return;
    else
        parentTempInfo=bioTree{parentTempInfo(1)}.node{parentTempInfo(2)}.parentNode;
    end
end
end
function divisionCount(bioTree,frameShift)
type3Node=0;
for iframe=frameShift+1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            if bioTree{iframe}.node{iNode}.type==3
                type3Node=type3Node+1;
            end
        end
    end
end
disp(strcat('type3Node=',num2str(type3Node)));
end


           
          