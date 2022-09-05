function bioTree=nodeEdgeReduction(bioTree,frameShift)
frameThreshold=15;
reducing=true;
i=1;
while reducing
[bioTree,nodeNumber1]=nodeShortNoiseFilter(bioTree,frameShift,frameThreshold) ;
[bioTree,edgeNumber1]=edgeRefine(bioTree,frameShift);
str=strcat(num2str(i),'.before Reduction: node Number =', num2str(nodeNumber1(1)),'; edge Number =',num2str(edgeNumber1(1)),'--after Redution: node Number =',num2str(nodeNumber1(2)),'; edge Number =',num2str(edgeNumber1(2)));
disp(str);
if nodeNumber1(1)==nodeNumber1(2) && edgeNumber1(1)==edgeNumber1(2)
    reducing=false;
end
i=i+1;
end
end