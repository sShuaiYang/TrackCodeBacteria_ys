% function newMatrix=directLinkNetworkForbinaryNetwork(linkMatrix,nodeInfo)
% linkMatrix=full(linkMatrix);
% newSize=size(linkMatrix,1);
% newMatrix=zeros(newSize);
% order=1:newSize;
% for i=1:size(linkMatrix,2)
%     iColumn=linkMatrix(:,i);
%     iOrder=order(iColumn);
%     for m=1:numel(iOrder)
%         for n=1:numel(iOrder)
%             if iOrder(m)~=iOrder(n) && nodeInfo(iOrder(n),5)~=nodeInfo(iOrder(m),5)
%                 newMatrix(iOrder(m),iOrder(n))=1;
%             end
%         end
%     end
% end
% end

function newMatrix=directLinkNetworkForbinaryNetwork(linkMatrix,nodeInfo)
linkMatrix=full(linkMatrix);
newSize=size(linkMatrix,1);
newMatrix=linkMatrix;
order=1:newSize;
for i=1:size(linkMatrix,2)
    for j=1:size(linkMatrix,2)
        if nodeInfo(i,5)==nodeInfo(j,5)
            newMatrix(i,j)=0;
        end
    end
end
end

function a 
degree=[];
distribution=[];
for i=1:size(newMatrix,1)
    iLine=newMatrix(i,:);
    degree(i)=numel(iLine(iLine~=0));
end
degreeMax=max(degree);
distribution(:,1)=1:degreeMax;
for i=1:degreeMax
distribution(i,2)=numel(degree(degree>=i));
end
end
function b
for i=1:numel(result{1})
result{1}(i).finalMatrix=directLinkNetworkForbinaryNetwork(result{1}(i).finalMatrix,clusterTree(i*200).nodeInfo);
end
for i=200:200:numel(bacteriaFrameInfo)
bacteriaInfo=bacteriaFrameInfo{i}.bacteriaInfo;
bacNum(i/200)=size(bacteriaInfo,1);
end
[bacteriaNum,expNum]=plotTransRateVsBacteriaNum(result);
end