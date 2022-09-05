function plotLinkTreeWithTwoChannel(linkMatrix,centroidInfo,leafYCoo,color,allList,bacteriaFrameInfo,bioTree)
% generally linkMatrix means the relationship for [N, L], N means all
% node, L means all leaf, centroidInfo here could be changed for the second
% column. The first column means the frame that the node begin\
% color=[1,0,1];
leafNum=numel(leafYCoo);
nodeNum=size(linkMatrix,1)-leafNum;
centroidInfo(nodeNum+1:end,2)=leafYCoo;
leafOrder=1:leafNum;
for i=1:nodeNum
    iLine=linkMatrix(nodeNum+1-i,nodeNum+1-i:end);
    centroid2=centroidInfo(:,2);
    centroid2=centroid2(nodeNum+1-i:end);
    leafPos=(centroid2(iLine~=0));
    centroidInfo(nodeNum+1-i,2)=mean(leafPos(~isnan(leafPos)));
end
% figure;
possibleNum=100:100:size(bioTree,2);
for i=1:size(linkMatrix,1)
    for j=i+1:size(linkMatrix,1)
        if linkMatrix(i,j)~=0
            linkTwoPointsAging(centroidInfo(i,:),centroidInfo(j,:),color,nodeNum,i,j,linkMatrix(i,j));
            if allList(j,3)==1
                nextNodeInfo=allList(j,:);
                isNode=bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{1}.isNode;
                if isNode==1
                    preNodeInfo=bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{1}.nodeInfo;
                    traceNum=size(bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.traceInfo.pixelIdxList,2);
                    for iframe=(possibleNum(possibleNum>=preNodeInfo(1) & possibleNum<=preNodeInfo(1)+traceNum-1))
                        bacInfo=[preNodeInfo(1),preNodeInfo(2),preNodeInfo(3),iframe-preNodeInfo(1)+1];
                        frameBacteriaInfo=bacteriaFrameInfo{iframe};
                        focusBac=frameBacteriaInfo.bacteriaInfo(:,1)==bacInfo(1) & frameBacteriaInfo.bacteriaInfo(:,2)==bacInfo(2) & frameBacteriaInfo.bacteriaInfo(:,3)==bacInfo(3) & frameBacteriaInfo.bacteriaInfo(:,4)==bacInfo(4);
                        rpf_gfpRatio=mean(getHigherData(frameBacteriaInfo.gfpImageInfo{focusBac}));
%                         rpf_gfpRatio=mean(frameBacteriaInfo.t_DimerImageInfo{focusBac}./frameBacteriaInfo.gfpImageInfo{focusBac});
                        text(iframe/1200,centroidInfo(j,2),num2str(roundn(rpf_gfpRatio,0)))
                    end
                end
                if isNode==0
                    preRootInfo=bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{1}.rootInfo;
                    traceNum=size(bioTree{preRootInfo(1)}.root{preRootInfo(2)}.traceInfo.pixelIdxList,2);
                    for iframe=(possibleNum(possibleNum>=preRootInfo(1) & possibleNum<=preRootInfo(1)+traceNum-1))
                        bacInfo=[preRootInfo(1),preRootInfo(2),0,iframe-preRootInfo(1)+1];
                        frameBacteriaInfo=bacteriaFrameInfo{iframe};
                        focusBac=frameBacteriaInfo.bacteriaInfo(:,1)==bacInfo(1) & frameBacteriaInfo.bacteriaInfo(:,2)==bacInfo(2) & frameBacteriaInfo.bacteriaInfo(:,3)==bacInfo(3) & frameBacteriaInfo.bacteriaInfo(:,4)==bacInfo(4);
                        rpf_gfpRatio=mean(getHigherData(frameBacteriaInfo.gfpImageInfo{focusBac}));
%                         rpf_gfpRatio=mean(frameBacteriaInfo.t_DimerImageInfo{focusBac}./frameBacteriaInfo.gfpImageInfo{focusBac});
                        text(iframe/1200,centroidInfo(j,2),num2str(roundn(rpf_gfpRatio,0)))
                    end
                end
            end
            if allList(j,3)==0
                nextLeafInfo=allList(j,:);
                is2Node=bioTree{nextLeafInfo(1)}.leavies{nextLeafInfo(2)}.is2Node;
                if is2Node==1
                    preNodeInfo=bioTree{nextLeafInfo(1)}.leavies{nextLeafInfo(2)}.nodeInfo;
                    traceNum=size(bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.traceInfo.pixelIdxList,2);
                    for iframe=(possibleNum(possibleNum>=preNodeInfo(1) & possibleNum<=preNodeInfo(1)+traceNum-1))
                        bacInfo=[preNodeInfo(1),preNodeInfo(2),preNodeInfo(3),iframe-preNodeInfo(1)+1];
                        frameBacteriaInfo=bacteriaFrameInfo{iframe};
                        focusBac=frameBacteriaInfo.bacteriaInfo(:,1)==bacInfo(1) & frameBacteriaInfo.bacteriaInfo(:,2)==bacInfo(2) & frameBacteriaInfo.bacteriaInfo(:,3)==bacInfo(3) & frameBacteriaInfo.bacteriaInfo(:,4)==bacInfo(4);
                        rpf_gfpRatio=mean(getHigherData(frameBacteriaInfo.gfpImageInfo{focusBac}));
%                         rpf_gfpRatio=sum(frameBacteriaInfo.t_DimerImageInfo{focusBac})/sum(frameBacteriaInfo.gfpImageInfo{focusBac});
%                         rpf_gfpRatio=mean(frameBacteriaInfo.t_DimerImageInfo{focusBac}./frameBacteriaInfo.gfpImageInfo{focusBac});

                        text(iframe/1200,centroidInfo(j,2),num2str(roundn(rpf_gfpRatio,0)))
                    end
                end
                if is2Node==0
                    preRootInfo=bioTree{nextLeafInfo(1)}.leavies{nextLeafInfo(2)}.rootInfo;
                    traceNum=size(bioTree{preRootInfo(1)}.root{preRootInfo(2)}.traceInfo.pixelIdxList,2);
                    for iframe=(possibleNum(possibleNum>=preRootInfo(1) & possibleNum<=preRootInfo(1)+traceNum-1))
                        bacInfo=[preRootInfo(1),preRootInfo(2),0,iframe-preRootInfo(1)+1];
                        frameBacteriaInfo=bacteriaFrameInfo{iframe};
                        focusBac=frameBacteriaInfo.bacteriaInfo(:,1)==bacInfo(1) & frameBacteriaInfo.bacteriaInfo(:,2)==bacInfo(2) & frameBacteriaInfo.bacteriaInfo(:,3)==bacInfo(3) & frameBacteriaInfo.bacteriaInfo(:,4)==bacInfo(4);
                        rpf_gfpRatio=mean(getHigherData(frameBacteriaInfo.gfpImageInfo{focusBac}));
%                         rpf_gfpRatio=sum(getHigherData(frameBacteriaInfo.t_DimerImageInfo{focusBac}))/sum(getHigherData(frameBacteriaInfo.gfpImageInfo{focusBac}));
%                         rpf_gfpRatio=mean(frameBacteriaInfo.t_DimerImageInfo{focusBac}./frameBacteriaInfo.gfpImageInfo{focusBac});
                        text(iframe/1200,centroidInfo(j,2),num2str(roundn(rpf_gfpRatio,0)))
                    end
                end
            end
        end
    end
end
end
function linkTwoPointsAging(c1,c2,color,nodeNum,i,j,para)
c1(1)=c1(1)/1200;
c2(1)=c2(1)/1200;
% support c1(x)<c2(x),link c1-c3-c2

% debug plus para(linkMatrix(i,j))
if para==-inf
    color=[0,0,0];
else
    colorAll=colormap(jet(8));
    color=colorAll(para,:);
end

if c1(1)>c2(1)
    c=c2;
    c2=c1;
    c1=c;
end
c3=[c1(1),c2(2)];
hold on;
line([c1(1),c3(1),c2(1)],[c1(2),c3(2),c2(2)],'Color',color,'LineWidth',1);
maker='.';
size=16;
if i==1
    plot(c1(1),c1(2),'Marker','^','MarkerSize',12,'Color',[1,0,0]);
else
    if j<=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
    end
    if j>=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        if c2(1)<12460/1200
            plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[85,107,47]/255);
        end
    end
end
end
function data=getHigherData(data)
data=sort(data);
number=numel(data);
data=data(end-ceil(0.8*number)+1:end);
end