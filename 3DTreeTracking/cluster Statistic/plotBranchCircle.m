function plotBranchCircle(linkMatrix)
branchList=linkMatrix.branchOrder(:,1);
branchOrder=linkMatrix.branchOrder(:,2);
branchListDiff=diff(branchList);
pointNum=size(branchList,1);
sita=pi*2/(pointNum-1);
colorNum=1;
figure;hold on;axis equal
for i=1:pointNum
    xdata=100*cos(sita*(i-1));
    ydata=100*sin(sita*(i-1));
    linkMatrix.branchOrder(i,3:4)=[xdata,ydata];
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
distMatrix=linkMatrix.distMatrix;
for i=1:size(distMatrix,1)
    for j=i+1:size(distMatrix,2)
        if distMatrix(i,j)==1
        line([linkMatrix.branchOrder(branchOrder(i),3),linkMatrix.branchOrder(branchOrder(j),3)],[linkMatrix.branchOrder(branchOrder(i),4),linkMatrix.branchOrder(branchOrder(j),4)]);hold on
        end
    end
end
end