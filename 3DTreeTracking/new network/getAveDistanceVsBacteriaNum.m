function [data] = getAveDistanceVsBacteriaNum( result )
for i=70:numel(result)
    clc
    disp(i)
    finalMatrix=result(i).finalMatrix;
    finalMatrix=full(finalMatrix);
    [data(i,2),data(i,1)]=getMaxCluster(finalMatrix);
end
end
function [averageDistance,num]=getMaxCluster(distMatrix)
gObj=biograph(distMatrix);
matrixSize=size(distMatrix,1);
[~,C]=conncomp(gObj,'directed','false');
c=sort(C);
cc=regionprops(c,'FilledArea');
maxArea=0;
for i=1:size(cc,1)
    if cc(i).FilledArea>maxArea
        maxArea=cc(i).FilledArea;
        orderNum=i;
    end
end
newMatrix=distMatrix(C==orderNum,C==orderNum);
if maxArea>=10
    gObj=biograph(newMatrix);
    newDist=allshortestpaths(gObj);
    averageDistance=mean(newDist(newDist~=inf & newDist~=0));
    num=size(newMatrix,1);
else
    averageDistance=+inf;
    num=0;
end
end

