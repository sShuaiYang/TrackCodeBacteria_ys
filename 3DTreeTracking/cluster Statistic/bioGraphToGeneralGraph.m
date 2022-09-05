function bioGraphToGeneralGraph(clusterTree,iframe,bacteriaFrameInfo)
[distMatrix,bacteriaList]=findLargestCLuster(full(clusterTree{iframe}.distMatrix),clusterTree{iframe}.bacteriaList);
% distMatrix=full(clusterTree{iframe}.distMatrix);
% bacteriaList=clusterTree{iframe}.bacteriaList;
bg=biograph(distMatrix);
% bg.LayoutType='hierarchical';
bg.edgeType='straight';
bg.LayoutType='equilibrium';
dolayout(bg);
for i=1:numel(bg.nodes)
    position(i,:)=bg.nodes(i).Position;
end
bacteriaList(:,6)=(bacteriaList(:,1)+bacteriaList(:,4))==iframe+1;
firstFrameBacteria=bacteriaFrameInfo{1}.bacteriaInfo(:,5);
firstFrameBacteria=unique(firstFrameBacteria);
bacteriaList(:,7)=ismember(bacteriaList(:,5),firstFrameBacteria);
bacteriaList(:,8)=sum(distMatrix,2);
maxLinkNum=max(bacteriaList(:,8));
minLinkNum=min(bacteriaList(:,8));
figure;hold on
for i=1:numel(bg.nodes)
    for j=i+1:numel(bg.nodes)
        if distMatrix(i,j)==1
            line([position(i,1),position(j,1)],[position(i,2),position(j,2)],'Color',[0.5,0.5,0.5]);
        end
    end
end
for i=1:numel(bg.nodes)
    iList=bacteriaList(i,:);
    if iList(7)==0
        filledColor='g';
    else
        filledColor='r';
    end
    if iList(6)==0
        filledColor=[1,1,0];
    end
    if iList(8)<=5
        objSize=7;
    else if iList(8)<=5+(maxLinkNum-minLinkNum)/3
            objSize=9;
        else if iList(8)<=5+(maxLinkNum-minLinkNum)/3*2
                objSize=11;
            else
                objSize=13;
            end
        end
    end
%     plot(position(i,1),position(i,2),'Marker',markerType,'MarkerFaceColor',filledColor,'MarkerSize',objSize)
plot(position(i,1),position(i,2),'MarkerFaceColor',filledColor,'MarkerSize',objSize,'Marker','o')
end
end

    