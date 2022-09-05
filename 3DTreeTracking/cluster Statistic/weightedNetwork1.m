% spatial size
function distMatrix=weightedNetwork1(clusterTree,bioTree,bacteriaFrameInfo,iframe)
lengthInfo=[];
for iFrame=1:size(bacteriaFrameInfo,2)
   lengthInfo=[lengthInfo;bacteriaFrameInfo{iFrame}.lengthInfo];
end
averageLength=mean(lengthInfo);
bacteriaList=clusterTree{iframe}.bacteriaList;
distMatrix=double(full(clusterTree{iframe}.distMatrix));
% properNode=(bacteriaList(:,1)+bacteriaList(:,4)-1)==numel(bioTree);
% distMatrix=distMatrix(properNode,properNode);
% bacteriaList=bacteriaList(properNode,:);
centroidInfo=zeros(size(distMatrix,1),2);
for i=1:size(distMatrix,1)
    if bacteriaList(i,3)==0
        centroidInfo(i,1:2)=bioTree{bacteriaList(i,1)}.root{bacteriaList(i,2)}.traceInfo.measurment{bacteriaList(i,4)}.Centroid;
    end
    if bacteriaList(i,3)~=0
        centroidInfo(i,1:2)=bioTree{bacteriaList(i,1)}.node{bacteriaList(i,2)}.Out{bacteriaList(i,3)}.traceInfo.measurment{bacteriaList(i,4)}.Centroid;
    end
end
distInfo=pdist2(centroidInfo,centroidInfo);
% radius=(bioTree{1}.imageSize(1).^2+bioTree{1}.imageSize(2).^2).^0.5/2;
% distInfo=distInfo/radius;
% distInfo(distInfo>=1)=1;
% distMatrix(distMatrix==1)=distInfo(distMatrix==1);
distInfo=distInfo/(averageLength*1.2);
distMatrix(distMatrix==1)=distInfo(distMatrix==1);
distMatrix=ceil(distMatrix);
end