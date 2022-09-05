function [largestMatrix,branchList]=findLargestCLuster(distMatrix,branchList)
gObj=biograph(distMatrix);
[~,C]=conncomp(gObj,'directed','false');
for i=1:max(C);
    clusterNum(i)=numel(C(C==i));
end
[~,iNum]=find(clusterNum==max(clusterNum));
largestMatrix=distMatrix(C==iNum,C==iNum);
branchList=branchList(C==iNum,:);
end
% distMatrix=full(distMatrix);
% BGobj = biograph(distMatrix);
% dolayout(BGobj);
% set(BGobj.nodes,'shape','circle','size',[80,80],'color',[1,0,0],'lineColor',[0,0,0])
% set(BGobj.edges,'LineColor',[0.5,0.5,0.5])