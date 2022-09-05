function testNodeTypeTracking(bioTree,nodeList)
for iList=1:size(nodeList,1)
    disp(iList);
    newTrace=fullNodeType1Tracking(bioTree,nodeList(iList,:));
end
end