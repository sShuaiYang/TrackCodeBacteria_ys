function tree=smallTreeStructureSeperate_cBF(bioTree)
bioTree1=bioTree;
branchList=bioTree{1}.branchList;
for i=1:size(branchList,1)
    linkInfo=[];
    bioTree=bioTreeCut(bioTree1,branchList(i,5));
    [linkMatrix,centroidInfo,leafList,leafNum,~,nodeInfo]=generateOneBranchPhytree(bioTreeCut(bioTree,branchList(i,5)),i,branchList(i,:),'normal');
    [linkRow,linkColumn]=find(linkMatrix==1);
    allList=[nodeInfo;leafList];
    for iLink=1:numel(linkRow)
        linkPre=allList(linkRow(iLink),:);
        if linkPre(4)==-1
            %             dataSize=numel(bioTree{linkPre(1)}.root{linkPre(2)}.CyOFP);
%             dataSize=numel(bioTree{linkPre(1)}.root{linkPre(2)}.linkTimer);
            linkInfo(iLink).pixelIdxList=bioTree{linkPre(1)}.root{linkPre(2)}.traceInfo.pixelIdxList;
            linkInfo(iLink).measurment=bioTree{linkPre(1)}.root{linkPre(2)}.traceInfo.measurment;
            %             linkInfo(iLink).CyOFP=bioTree{linkPre(1)}.root{linkPre(2)}.CyOFP;
            %             linkInfo(iLink).GFP=bioTree{linkPre(1)}.root{linkPre(2)}.GFP;
            %             linkInfo(iLink).RFP=bioTree{linkPre(1)}.root{linkPre(2)}.RFP;
            %             linkInfo(iLink).mScalet=bioTree{linkPre(1)}.root{linkPre(2)}.mScalet;
            %             linkInfo(iLink).CyPet=bioTree{linkPre(1)}.root{linkPre(2)}.CyPet;
            %             linkInfo(iLink).Venus=bioTree{linkPre(1)}.root{linkPre(2)}.Venus;
            %             linkInfo(iLink).mAmetrine=bioTree{linkPre(1)}.root{linkPre(2)}.mAmetrine;
            linkInfo(iLink).linkTimer=bioTree{linkPre(1)}.root{linkPre(2)}.linkTimer;
            %             linkInfo(iLink).linkAge=ones(1,numel(bioTree{linkPre(1)}.root{linkPre(2)}.CyOFP));
            %             linkInfo(iLink).redControl=bioTree{linkPre(1)}.root{linkPre(2)}.RedControl;
            %             linkInfo(iLink).blueControl=bioTree{linkPre(1)}.root{linkPre(2)}.BlueControl;
            %             linkInfo(iLink).greenControl=bioTree{linkPre(1)}.root{linkPre(2)}.GreenControl;
%             length=ones(dataSize,1)*nan;
%             growthRate=ones(dataSize,1)*nan;
            for iLine=2:2:numel(bioTree{linkPre(1)}.root{linkPre(2)}.traceInfo.measurment)
                length(iLine)=bioTree{linkPre(1)}.root{linkPre(2)}.traceInfo.measurment{iLine}.MajorAxisLength;
            end
%             for iLine=2:2:numel(bioTree{linkPre(1)}.root{linkPre(2)}.traceInfo.measurment)-2
%                 growthRate(iLine)=log2(length(iLine+2)/length(iLine))/(linkInfo(iLink).linkTimer(iLine+2)-linkInfo(iLink).linkTimer(iLine));
%             end
%             if isempty(iLine)
%                 growthRate(2)=nan;
%             else
%                 growthRate(iLine+2)=growthRate(iLine);
%             end
%             linkInfo(iLink).growthRate=growthRate;
            %             linkInfo(iLink).fpProduceCyOFP=getProduceRate(iLine,linkInfo(iLink).CyOFP,growthRate,linkInfo(iLink).linkTimer,dataSize);
            %             linkInfo(iLink).fpProduceGFP=getProduceRate(iLine,linkInfo(iLink).GFP,growthRate,linkInfo(iLink).linkTimer,dataSize);
            %             linkInfo(iLink).fpProduceRFP=getProduceRate(iLine,linkInfo(iLink).RFP,growthRate,linkInfo(iLink).linkTimer,dataSize);
            %             linkInfo(iLink).fpProducemScalet=getProduceRate(iLine,linkInfo(iLink).mScalet,growthRate,linkInfo(iLink).linkTimer,dataSize);
            %             linkInfo(iLink).fpProduceCyPet=getProduceRate(iLine,linkInfo(iLink).CyPet,growthRate,linkInfo(iLink).linkTimer,dataSize);
            %             linkInfo(iLink).fpProduceVenus=getProduceRate(iLine,linkInfo(iLink).Venus,growthRate,linkInfo(iLink).linkTimer,dataSize);
            %             linkInfo(iLink).fpProducemAmetrine=getProduceRate(iLine,linkInfo(iLink).mAmetrine,growthRate,linkInfo(iLink).linkTimer,dataSize);
        end
        if linkPre(4)==1
            nextInfo=[];
            linkNext=allList(linkColumn(iLink),:);
            for iOut=1:numel(bioTree{linkPre(1)}.node{linkPre(2)}.Out)
                if bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.is2Node==1
                    nextInfo=[bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.nodeInfo,1];
                    if isequal(linkNext,nextInfo)
                        break
                    end
                else
                    nextInfo=[bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.leafInfo,0,0];
                    if isequal(linkNext,nextInfo)
                        break
                    end
                end
            end
            %             dataSize=numel(bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.CyOFP);
%             dataSize=numel(bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.linkTimer);
            linkInfo(iLink).pixelIdxList=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.traceInfo.pixelIdxList;
            linkInfo(iLink).measurment=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.traceInfo.measurment;
            %             linkInfo(iLink).CyOFP=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.CyOFP;
            %             linkInfo(iLink).GFP=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.GFP;
            %             linkInfo(iLink).RFP=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.RFP;
            %             linkInfo(iLink).mScalet=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.mScalet;
            %             linkInfo(iLink).CyPet=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.CyPet;
            %             linkInfo(iLink).Venus=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.Venus;
            %             linkInfo(iLink).mAmetrine=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.mAmetrine;
            %             linkInfo(iLink).linkTimer=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.linkTimer;
            %             linkInfo(iLink).linkAge=ones(1,numel(bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.CyOFP))*max(bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.isOld1,bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.isOld2);
            %             linkInfo(iLink).redControl=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.RedControl;
            %             linkInfo(iLink).blueControl=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.BlueControl;
            %             linkInfo(iLink).greenControl=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.GreenControl;
%             length=ones(dataSize,1)*nan;
%             growthRate=ones(dataSize,1)*nan;
            for iLine=2:2:numel(bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.traceInfo.measurment)
                length(iLine)=bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.traceInfo.measurment{iLine}.MajorAxisLength;
            end
%             for iLine=2:2:numel(bioTree{linkPre(1)}.node{linkPre(2)}.Out{iOut}.traceInfo.measurment)-2
%                 growthRate(iLine)=log2(length(iLine+2)/length(iLine))/(linkInfo(iLink).linkTimer(iLine+2)-linkInfo(iLink).linkTimer(iLine));
%             end
%             if isempty(iLine)
%                 growthRate(2)=0;
%             else
%                 growthRate(iLine+2)=growthRate(iLine);
%             end
%             linkInfo(iLink).growthRate=growthRate;
            %             linkInfo(iLink).fpProduceCyOFP=getProduceRate(iLine,linkInfo(iLink).CyOFP,growthRate,linkInfo(iLink).linkTimer,dataSize);
            %             linkInfo(iLink).fpProduceGFP=getProduceRate(iLine,linkInfo(iLink).GFP,growthRate,linkInfo(iLink).linkTimer,dataSize);
            %             linkInfo(iLink).fpProduceRFP=getProduceRate(iLine,linkInfo(iLink).RFP,growthRate,linkInfo(iLink).linkTimer,dataSize);
            %             linkInfo(iLink).fpProducemScalet=getProduceRate(iLine,linkInfo(iLink).mScalet,growthRate,linkInfo(iLink).linkTimer,dataSize);
            %             linkInfo(iLink).fpProduceCyPet=getProduceRate(iLine,linkInfo(iLink).CyPet,growthRate,linkInfo(iLink).linkTimer,dataSize);
            %             linkInfo(iLink).fpProduceVenus=getProduceRate(iLine,linkInfo(iLink).Venus,growthRate,linkInfo(iLink).linkTimer,dataSize);
            %             linkInfo(iLink).fpProducemAmetrine=getProduceRate(iLine,linkInfo(iLink).mAmetrine,growthRate,linkInfo(iLink).linkTimer,dataSize);
        end
    end
    tree(i).linkInfo=linkInfo;
    tree(i).linkRow=linkRow;
    tree(i).linkColumn=linkColumn;
    tree(i).linkMatrix=linkMatrix;
    tree(i).centroidInfo=centroidInfo;
    tree(i).leafNum=leafNum;
    tree(i).cutNum=bioTree{1}.bioTreeTimer(branchList(i,5));
    tree(i).allList=allList;
    tree(i).timer=bioTree{1}.bioTreeTimer(allList(:,1))';
    %     tree(i).tag=bioTree{1}.mip.feildTag.tag{bioTree{1}.fieldNum};
    %     tree(i).tagValue=bioTree{1}.mip.feildTag.tagValue(bioTree{1}.fieldNum);
    %     tree(i).scale=bioTree{1}.mip.Calib.scale;
    tree(i).fieldNum=bioTree{1}.fieldNum;
end
end
function fpProduce=getProduceRate(iLine,fluoFP,growthRate,linkTimer,dataSize)
fpProduce=ones(dataSize,1)*nan;
if isempty(iLine)
    fpProduce(2)=nan;
else
    gfpDataIndex=find(~isnan(fluoFP));
    if numel(gfpDataIndex)~=1
        for iGFP=1:numel(gfpDataIndex)-1
            fpProduce(gfpDataIndex(iGFP))=growthRate(gfpDataIndex(iGFP))*fluoFP(gfpDataIndex(iGFP))+(fluoFP(gfpDataIndex(iGFP+1))-fluoFP(gfpDataIndex(iGFP)))/(linkTimer(gfpDataIndex(iGFP+1))-linkTimer(gfpDataIndex(iGFP)));
        end
        fpProduce(gfpDataIndex(iGFP+1))=fpProduce(gfpDataIndex(iGFP));
    end
end
end