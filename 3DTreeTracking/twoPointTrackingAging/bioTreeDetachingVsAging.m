function ageFate=bioTreeDetachingVsAging(bioTree)
ageFate=[];
limitRegionIdx=getLimitRegion(bioTree{1}.imageProcessingInfo.cropInfo);
for iframe=1:numel(bioTree)
    for iNode=1:size(bioTree{iframe}.node,2)
        age=[];
        fate=[];
        n=0;
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            if bioTree{iframe}.node{iNode}.Out{iOut}.isOld1>=0 && bioTree{iframe}.node{iNode}.Out{iOut}.isOld2>=0
                n=n+1;
                age(n)=max(bioTree{iframe}.node{iNode}.Out{iOut}.isOld1,bioTree{iframe}.node{iNode}.Out{iOut}.isOld2);
                fate(n)=bioTree{iframe}.node{iNode}.Out{iOut}.is2Node;
                if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==0
                    if bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo(1)==numel(bioTree)
                        age(n)=[];
                        fate(n)=[];
                        n=n-1;
                        continue
                    end
                    leafInfo=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
                    leafPixelDetail=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leafPixelDetail;
                    if any(ismember(leafPixelDetail,limitRegionIdx))
                        age(n)=[];
                        fate(n)=[];
                        n=n-1;
                        continue
                    end
                end
            end
        end
        if ~(numel(age)==2 && age(1)==2 && age(2)==2)
            for i=1:numel(age)
                if age(i)==1
                    ageFate=[ageFate;0,fate(i)];
                end
                if age(i)>1
                    ageFate=[ageFate;1,fate(i)];
                end
            end
        end
    end
end
end
function limitRegionIdx=getLimitRegion(cropInfo)
cropInfo1=bwmorph(cropInfo,'remove');
cropInfo1=imdilate(cropInfo1,true(13));
limitRegion=cropInfo1 & cropInfo;
cc=bwconncomp(limitRegion);
limitRegionIdx=cc.PixelIdxList{1,1};
end
function testtest
u=[];
for i=1:2
age=ageFate(:,1);
fate=ageFate(:,2);
ageNum=numel(age(age==(i-1)));
fateAge=fate(age==(i-1));
detachFate=numel(fateAge(fateAge==0));
u=[u;ageNum,detachFate];
end;ratio=u(:,2)./u(:,1)
end
                