function newinfo=biomatrixSwallow(bacteriaFrameInfo,clusterTree,Swallowprobability)
bacteriaFrameInfo=bacteriaFrameInfo(201:200:end);
clusterTree=clusterTree(201:200:end);
InfoSize=size(bacteriaFrameInfo,2);
for iNum=1:InfoSize
    biotreeInfo{iNum}.frameNum=iNum*200+1;
    biotreeInfo{iNum}.bacteriaInfo=bacteriaFrameInfo{iNum}.bacteriaInfo;
    biotreeInfo{iNum}.allbacteriaList=clusterTree{iNum}.bacteriaList;
    biotreeInfo{iNum}.bacnum=size(biotreeInfo{iNum}.bacteriaInfo,1);
    biotreeInfo{iNum}.allbacnum=size(biotreeInfo{iNum}.allbacteriaList,1);
    biotreeInfo{iNum}.centroidInfo=bacteriaFrameInfo{iNum}.centroidInfo;
    biotreeInfo{iNum}.avelength=mean(bacteriaFrameInfo{iNum}.lengthInfo);
    biotreeInfo{iNum}.biomatrix=clusterTree{iNum}.distMatrix;
    comparelist=ones(biotreeInfo{iNum}.allbacnum,1);
    for ibac=1:biotreeInfo{iNum}.bacnum
        num=biotreeInfo{iNum}.bacteriaInfo(ibac,1)+1e5*biotreeInfo{iNum}.bacteriaInfo(ibac,2)+1e7*biotreeInfo{iNum}.bacteriaInfo(ibac,3);
        infonum=find(biotreeInfo{iNum}.allbacteriaList(1:biotreeInfo{iNum}.allbacnum,1)+1e5*biotreeInfo{iNum}.allbacteriaList(1:biotreeInfo{iNum}.allbacnum,2)+1e7*biotreeInfo{iNum}.allbacteriaList(1:biotreeInfo{iNum}.allbacnum,3)==num*comparelist);
        biotreeInfo{iNum}.bacteriaInfo(ibac,6)=infonum;
        biotreeInfo{iNum}.allbacteriaList(infonum,6)=ibac;
    end
end
for iNum=1:InfoSize-1
    disp(iNum);
    comparelist=ones(biotreeInfo{iNum+1}.allbacnum,1);
    for ibac=1:biotreeInfo{iNum}.allbacnum
        if biotreeInfo{iNum}.allbacteriaList(ibac,1)~=biotreeInfo{iNum+1}.allbacteriaList(ibac,1)
            num=biotreeInfo{iNum+1}.allbacteriaList(ibac,1)+1e5*biotreeInfo{iNum+1}.allbacteriaList(ibac,2)+1e7*biotreeInfo{iNum+1}.allbacteriaList(ibac,3);
            infonum=find(biotreeInfo{iNum+1}.allbacteriaList(1:biotreeInfo{iNum+1}.allbacnum,1)+1e5*biotreeInfo{iNum+1}.allbacteriaList(1:biotreeInfo{iNum+1}.allbacnum,2)+1e7*biotreeInfo{iNum+1}.allbacteriaList(1:biotreeInfo{iNum+1}.allbacnum,3)==num*comparelist+1e7);
            if ~isempty(infonum)
                biotreeInfo{iNum}.allbacteriaList(ibac,7)=infonum;
            else
                biotreeInfo{iNum}.allbacteriaList(ibac,7)=0;
            end
        else
            biotreeInfo{iNum}.allbacteriaList(ibac,7)=0;
        end
    end
    biotreeInfo{iNum}.allbacteriaList(1:biotreeInfo{iNum}.allbacnum,8)=0;
    biotreeInfo{iNum}.bacteriaInfo(1:biotreeInfo{iNum}.bacnum,8)=0;
end
for iNum=1:InfoSize
    disp(iNum);
    for ibac=1:biotreeInfo{iNum}.allbacnum
        if biotreeInfo{iNum}.allbacteriaList(ibac,8)==1
            if biotreeInfo{iNum}.allbacteriaList(ibac,6)~=0
               biotreeInfo{iNum}.bacteriaInfo(biotreeInfo{iNum}.allbacteriaList(ibac,6),8)=1;
            end
        end
    end
    biotreeInfo{iNum}.swallowList=[];
    properclusterdis2=(1.2*biotreeInfo{iNum}.avelength)^2;
    distance2List=ones(biotreeInfo{iNum}.bacnum,1).*properclusterdis2;
    for ibac=1:biotreeInfo{iNum}.bacnum
        deltaposition=biotreeInfo{iNum}.centroidInfo;
        for i=1:biotreeInfo{iNum}.bacnum
            deltaposition(i,1)=deltaposition(i,1)-biotreeInfo{iNum}.centroidInfo(ibac,1);
            deltaposition(i,2)=deltaposition(i,2)-biotreeInfo{iNum}.centroidInfo(ibac,2);
            if biotreeInfo{iNum}.bacteriaInfo(ibac,8)==1
                deltaposition(i,1:2)=[1000 1000];
            end
            biotreeInfo{iNum}.bacteriaInfo(ibac,7)=size(find(deltaposition(1:biotreeInfo{iNum}.bacnum,1).^2+deltaposition(1:biotreeInfo{iNum}.bacnum,2).^2<=distance2List),1);
        end
        if biotreeInfo{iNum}.bacteriaInfo(ibac,8)~=1
            isswallowed=swallowjudge(biotreeInfo{iNum}.bacteriaInfo(ibac,7),Swallowprobability);
             if isswallowed
                 biotreeInfo{iNum}.bacteriaInfo(ibac,8)=1;
                 biotreeInfo{iNum}.allbacteriaList(biotreeInfo{iNum}.bacteriaInfo(ibac,6),8)=1;
             end
        end        
    end
    for ibac=1:biotreeInfo{iNum}.allbacnum 
        if biotreeInfo{iNum}.allbacteriaList(ibac,8)==1 
            biotreeInfo{iNum}.swallowList(end+1)=ibac;
            if iNum==InfoSize
                continue
            end
            biotreeInfo{iNum+1}.allbacteriaList(ibac,8)=1;
            if biotreeInfo{iNum}.allbacteriaList(ibac,7)~=0
                biotreeInfo{iNum+1}.allbacteriaList(biotreeInfo{iNum}.allbacteriaList(ibac,7),8)=1;
            end
        end
    end
    biotreeInfo{iNum}.newmatrix=biotreeInfo{iNum}.biomatrix;
    if ~isempty(biotreeInfo{iNum}.swallowList)
        for iswa=size(biotreeInfo{iNum}.swallowList,2):-1:1
        biotreeInfo{iNum}.newmatrix=deleteaLineinMatrix(biotreeInfo{iNum}.newmatrix,biotreeInfo{iNum}.swallowList(iswa));
        end
    end
end
newinfo=biotreeInfo;
end
function isswallowed=swallowjudge(clustersize,Swallowprobability)
Swallowprobability=Swallowprobability/clustersize;
if rand(1)<=Swallowprobability
    isswallowed=1;
else
    isswallowed=0;
end
end
function newma=deleteaLineinMatrix(oldma,linenum)
newma=oldma;
newma(:,linenum)=[];
newma(linenum,:)=[];
end

