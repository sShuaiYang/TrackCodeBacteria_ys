function adjacencyMatrix=generateRandomGraph(nodeNum,edgeNum)
adjacencyMatrix=false(nodeNum);
linkIdx=getRandomLink(edgeNum,nodeNum);
XYIdx=recTangleIdxToXYIndex(linkIdx);
for i=1:size(XYIdx,1)
    adjacencyMatrix(XYIdx(i,1),XYIdx(i,2))=true;
    adjacencyMatrix(XYIdx(i,2),XYIdx(i,1))=true;
end
end
function XYIdx=recTangleIdxToXYIndex(recTangleIdx)
nSeries=fix((1+(1+8*recTangleIdx).^0.5)/2);
ySeries=recTangleIdx-nSeries.*(nSeries-1)/2;
ySeries(ySeries==0)=nSeries(ySeries==0);
XYIdx=cat(2,nSeries+1,ySeries);
end
function linkIdx=getRandomLink(edgeNum,nodeNum)
targetNum=1:nodeNum*(nodeNum-1)/2;
linkIdx=[];
for i=1:edgeNum
    randomNum=fix(rand*numel(targetNum));
    linkIdx=[linkIdx;targetNum(randomNum)];
    targetNum(randomNum)=[];
end
end
    