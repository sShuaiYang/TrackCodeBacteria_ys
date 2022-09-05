% iBac=9;
% imageStack1=imageStack;
% ibacteria=u.bacteriaList(iBac,:);
% if ibacteria(3)==0
% centroid1=bioTree{ibacteria(1)}.root{ibacteria(2)}.traceInfo.measurment{ibacteria{4}}.Centroid;
% imageStack1(bioTree{ibacteria(1)}.root{ibacteria(2)}.traceInfo.pixelIdxList{ibacteria(4)})=255;
% end
% if ibacteria(3)==1
% centroid1=bioTree{ibacteria(1)}.node{ibacteria(2)}.Out{ibacteria(3)}.traceInfo.measurment{ibacteria(4)}.Centroid;
% imageStack1(bioTree{ibacteria(1)}.node{ibacteria(2)}.Out{ibacteria(3)}.traceInfo.pixelIdxList{ibacteria(4)})=255;
% end
% imshow(imageStack1)
% hold on;plot(centroid1(1),centroid1(2),'yo')
% for i=1:2660
% distMatrix=full(u.distMatrix);
% iLine=distMatrix(i,:);
% degree=numel(iLine(iLine==1));
% ibacteria=u.bacteriaList(i,:);
% if ibacteria(3)==0
% centroid=bioTree{ibacteria(1)}.root{ibacteria(2)}.traceInfo.measurment{ibacteria{4}}.Centroid;
% end
% if ibacteria(3)==1
% centroid=bioTree{ibacteria(1)}.node{ibacteria(2)}.Out{ibacteria(3)}.traceInfo.measurment{ibacteria(4)}.Centroid;
% end
% text(centroid(1),centroid(2),num2str(degree))
% end

% test clusterTree & distMatrix
for i=1:size(bacteriaInfo,1)
    ibac=bacteriaInfo(i,:);
    if ibac(1)+ibac(4)-1~=12460
        if ibac(3)==0
            is2Node=bioTree{ibac(1)}.root{ibac(2)}.is2Node;
        else
            is2Node=bioTree{ibac(1)}.node{ibac(2)}.Out{ibac(3)}.is2Node;
        end
    end
    if is2Node~=0
        disp(ibac)
    end
end

% link high degree bacteria
aimI=1120;
aimIbac=clusterTree{end}.bacteriaList(aimI,:);
if aimIbac(3)==0
    centroidAim=bioTree{aimIbac(1)}.root{aimIbac(2)}.traceInfo.measurment{aimIbac(4)}.Centroid;
else
    centroidAim=bioTree{aimIbac(1)}.node{aimIbac(2)}.Out{aimIbac(3)}.traceInfo.measurment{aimIbac(4)}.Centroid;
end
iLine=clusterTree{end}.distMatrix(aimI,:);
hold on
for iii=1:numel(iLine)
    if iLine(iii)==1
        linkBac=clusterTree{end}.bacteriaList(iii,:);
        if linkBac(3)==0
            centroidLink=bioTree{linkBac(1)}.root{linkBac(2)}.traceInfo.measurment{linkBac(4)}.Centroid;
        else
            centroidLink=bioTree{linkBac(1)}.node{linkBac(2)}.Out{linkBac(3)}.traceInfo.measurment{linkBac(4)}.Centroid;
        end
        line([centroidAim(1),centroidLink(1)],[centroidAim(2),centroidLink(2)],'Color',[1,1,0]);
    end
end