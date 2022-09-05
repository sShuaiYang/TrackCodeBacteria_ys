bacteriaList=clusterTree{end}.bacteriaList;
distMatrix=clusterTree{end}.distMatrix;
detachingOne=(bacteriaList(:,1)+bacteriaList(:,4)-1)~=numel(bacteriaFrameInfo);
newMatrix=distMatrix(detachingOne,detachingOne);
degree=zeros(size(newMatrix,1),1);
for i=1:size(newMatrix,1)
    iLine=newMatrix(i,:);
    degree(i)=numel(iLine(iLine~=0));
end
degree(degree==0)=[];
maxDegree=max(degree);
regionHist(1,:)=1:maxDegree;
for i=1:maxDegree
    regionHist(2,i)=numel(degree(degree>=i));
end