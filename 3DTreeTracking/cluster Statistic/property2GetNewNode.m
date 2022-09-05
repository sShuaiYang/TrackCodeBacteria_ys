function property2GetNewNode(matrixA,matrixB)
getNodeNum=[];
for i=1:size(matrixA,1)
    iALine=matrixA(i,:);
    iBLine=matrixB(i,:);
    getNodeNum=[getNodeNum;[numel(iALine(iALine==1)),numel(iBLine(iBLine==1))-numel(iALine(iALine==1))]];
end
getNodeNum(getNodeNum(:,1)==0,:)=[];
[~,sortDegree]=sort(getNodeNum(:,1));
getNodeNum=getNodeNum(sortDegree,:);
i=2;
while ~(i==size(getNodeNum,1) && getNodeNum(i,1)~=getNodeNum(i-1,1))
    if getNodeNum(i,1)==getNodeNum(i-1,1)
        getNodeNum(i-1,2)=getNodeNum(i-1,2)+getNodeNum(i,2);
        getNodeNum(i,:)=[];
        i=i-1;
    end
    i=i+1;
end
getNodeNum(:,2)=getNodeNum(:,2)/sum(getNodeNum(:,2));
accumulateDis=getNodeNum;
for i=2:size(getNodeNum,1)
    accumulateDis(i,2)=sum(getNodeNum(1:i,2));
end
end