function plotCluster_branchCircle(NBCMatrix)
branchList=NBCMatrix(:,2);
branchListDiff=diff(branchList);
% branchNum=numel(branchListDiff(branchListDiff==1))+1;
% myColor=colormap(jet(branchNum));
pointNum=size(branchList,1);
sita=pi*2/(pointNum-1);
colorNum=1;
figure;hold on;axis equal
for i=1:pointNum
    xdata=100*cos(sita*(i-1));
    ydata=100*sin(sita*(i-1));
    NBCMatrix(i,4:5)=[xdata,ydata];
    if i>1 && branchListDiff(i-1)==1
        colorNum=colorNum+1;
    end
    if mod(colorNum,2)==0
    plot(xdata,ydata,'Color',[1,0,0],'Marker','o','LineStyle','none');hold on
    end
    if mod(colorNum,2)==1
    plot(xdata,ydata,'Color',[0,1,0],'Marker','o','LineStyle','none');hold on
    end
end
clusterList=NBCMatrix(:,3);
[~,clusterList]=sort(clusterList);
NBCMatrix=NBCMatrix(clusterList,:);
for iCluster=1:NBCMatrix(end,3)
    properOne=NBCMatrix(NBCMatrix(:,3)==iCluster,4:5);
    for i=1:size(properOne,1)
        for j=i+1:size(properOne,1)
            line([properOne(i,1),properOne(j,1)],[properOne(i,2),properOne(j,2)]);hold on
        end
    end
end
end
