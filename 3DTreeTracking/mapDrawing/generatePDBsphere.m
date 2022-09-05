% function [treeModel,centroidInfoAll]=generatePDBsphere(oriStructure,bacteriaFrameInfo)
function [treeModel]=generatePDBsphere(oriStructure,bacteriaFrameInfo,bioTree)
orderBranchList=[6];
for iframe=1:5:numel(bacteriaFrameInfo)
    branchList=bacteriaFrameInfo{iframe}.bacteriaInfo(:,5);
    needBranch=ismember(branchList,orderBranchList);
    bacteriaFrameInfo{iframe}.centroidInfo=bacteriaFrameInfo{iframe}.centroidInfo(needBranch,:);
    bacteriaFrameInfo{iframe}.bacteriaInfo=bacteriaFrameInfo{iframe}.bacteriaInfo(needBranch,:);
end
element{1}='Li';
element{2}='He';
element{3}='Be';
element{4}='B';
element{5}='C';
element{6}='N';
element{7}='O';
element{8}='F';
element{9}='Ne';
element{10}='Na';
centroidInfoAll=[];
n=0;
xyDiv=3;
zDiv=40;
for iBranch=1:numel(orderBranchList)
    for iframe=1:10:numel(bacteriaFrameInfo)
        %     for iframe=2000:5:4500
        bacteriaInfo=bacteriaFrameInfo{iframe}.bacteriaInfo(:,5)==orderBranchList(iBranch);
        centroidInfo=bacteriaFrameInfo{iframe}.centroidInfo(bacteriaInfo,:);
%         centroidInfo=centroidInfo(centroidInfo(:,1)>=1800,:);
%         centroidInfo=centroidInfo(centroidInfo(:,2)<=1500 & centroidInfo(:,2)>=650,:);
        centroidInfoAll=[centroidInfoAll;centroidInfo];
        for iCen=1:size(centroidInfo,1)
            n=n+1;
            treeModel.Model.Atom(n)=oriStructure;
            treeModel.Model.Atom(n).X=centroidInfo(iCen,1)/xyDiv;
            treeModel.Model.Atom(n).Y=centroidInfo(iCen,2)/xyDiv;
            treeModel.Model.Atom(n).Z=double(iframe)/zDiv;
            treeModel.Model.Atom(n).element=element{iBranch};
            treeModel.Model.Atom(n).AtomNameStruct.occupancy=1;
            treeModel.Model.Atom(n).AtomSerNo=n;
        end
    end
end
% add differentNode;
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
[bioTree,branchList,allList]=divisionFinder(bioTree,branchList);
for iBranch=1:numel(orderBranchList)
    focusBranch=orderBranchList(iBranch);
    aimNode=branchList(focusBranch,:);
    allRoot=bioTree{aimNode(1)}.node{aimNode(2)}.allRoot;
    centroidInfo=bioTree{allRoot(1)}.root{allRoot(2)}.rootMeasurment.Centroid;
%     centroidInfo=centroidInfo(centroidInfo(:,1)>=1800,:);
%     centroidInfo=centroidInfo(centroidInfo(:,2)<=1500 & centroidInfo(:,2)>=650,:);
    if ~isempty(centroidInfo)
        n=n+1;
        treeModel.Model.Atom(n)=oriStructure;
        treeModel.Model.Atom(n).X=centroidInfo(1,1)/xyDiv;
        treeModel.Model.Atom(n).Y=centroidInfo(1,2)/xyDiv;
        treeModel.Model.Atom(n).Z=double(allRoot(1))/zDiv;
        treeModel.Model.Atom(n).element='Cu';
        treeModel.Model.Atom(n).AtomNameStruct.occupancy=4;
        treeModel.Model.Atom(n).AtomSerNo=n;
    end
    allLeaf=bioTree{aimNode(1)}.node{aimNode(2)}.allLeaf;
    for iLeaf=1:size(allLeaf,1)
        eachLeaf=allLeaf(iLeaf,:);
        if eachLeaf(1)~=numel(bioTree)
            centroidInfo=bioTree{eachLeaf(1)}.leavies{eachLeaf(2)}.leafMeasurment.Centroid;
%             centroidInfo=centroidInfo(centroidInfo(:,1)>=1800,:);
%             centroidInfo=centroidInfo(centroidInfo(:,2)<=1500 & centroidInfo(:,2)>=650,:);
            if ~isempty(centroidInfo)
                n=n+1;
                treeModel.Model.Atom(n)=oriStructure;
                treeModel.Model.Atom(n).X=centroidInfo(1,1)/xyDiv;
                treeModel.Model.Atom(n).Y=centroidInfo(1,2)/xyDiv;
                treeModel.Model.Atom(n).Z=double(eachLeaf(1))/zDiv;
                treeModel.Model.Atom(n).element='Hg';
                treeModel.Model.Atom(n).AtomNameStruct.occupancy=4;
                treeModel.Model.Atom(n).AtomSerNo=n;
            end
        end
    end
    allNode=bioTree{aimNode(1)}.node{aimNode(2)}.allNode;
    for iNode=1:size(allNode,1)
        eachNode=allNode(iNode,:);
        if bioTree{eachNode(1)}.node{eachNode(2)}.In{1}.isNode==1
            preNode=bioTree{eachNode(1)}.node{eachNode(2)}.In{1}.nodeInfo;
            centroidInfo=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.measurment{end}.Centroid;
        end
        if bioTree{eachNode(1)}.node{eachNode(2)}.In{1}.isNode==0
            preRoot=bioTree{eachNode(1)}.node{eachNode(2)}.In{1}.rootInfo;
            centroidInfo=bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.measurment{end}.Centroid;
        end
        %         centroidInfo=centroidInfo(centroidInfo(:,1)>=800,:);
        %         centroidInfo=centroidInfo(centroidInfo(:,2)<=1600,:);
%         centroidInfo=centroidInfo(centroidInfo(:,1)>=1800,:);
%         centroidInfo=centroidInfo(centroidInfo(:,2)<=1500 & centroidInfo(:,2)>=650,:);
        if ~isempty(centroidInfo)
            n=n+1;
            treeModel.Model.Atom(n)=oriStructure;
            treeModel.Model.Atom(n).X=centroidInfo(1,1)/xyDiv;
            treeModel.Model.Atom(n).Y=centroidInfo(1,2)/xyDiv;
            treeModel.Model.Atom(n).Z=double(eachNode(1)-1)/zDiv;
            treeModel.Model.Atom(n).element='Ag';
            treeModel.Model.Atom(n).AtomNameStruct.occupancy=4;
            treeModel.Model.Atom(n).AtomSerNo=n;
        end
    end
end
end
