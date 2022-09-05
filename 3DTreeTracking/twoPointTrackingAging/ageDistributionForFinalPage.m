function ageList=ageDistributionForFinalPage(bioTree,bacteriaFrameInfo)
finalInfo=bacteriaFrameInfo{end};
ageList=[];
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
for i=1:size(finalInfo.bacteriaInfo,1)
    bacInfo=finalInfo.bacteriaInfo(i,:);
    branchIndex=bacInfo(5);
    branchInfo=branchList(branchIndex,:);
    if branchInfo(3)~=0
        allRoot=bioTree{branchInfo(1)}.node{branchInfo(2)}.allRoot;
        if numel(bioTree)-allRoot(1)>=1200*8
            if bacInfo(3)==0
                p1Old=bioTree{bacInfo(1)}.root{bacInfo(2)}.isOld1;
                p2Old=bioTree{bacInfo(1)}.root{bacInfo(2)}.isOld2;
            end
            if bacInfo(3)~=0
                p1Old=bioTree{bacInfo(1)}.node{bacInfo(2)}.Out{bacInfo(3)}.isOld1;
                p2Old=bioTree{bacInfo(1)}.node{bacInfo(2)}.Out{bacInfo(3)}.isOld2;
            end
            if p1Old~=-inf && p2Old~=-inf;
                ageList=[ageList;max(p1Old,p2Old)];
            end
        end
    end
end
end
    