function GrowthRateVSratioInDiffAveGrowthRate
% you must get the variaty aveGrowthRateVSgrowthRate
binSpaceX1=logspace(-2,-1,21);
binSpaceY1=logspace(-6,0,61);
[countMap1,~,~]=get2Dhist(aveGrowthRateVSgrowthRate(:,1),aveGrowthRateVSgrowthRate(:,2),binSpaceX1,binSpaceY1);
yBin1=binSpaceY1(2:2:end);
for i=1:30
    countMap1(:,i)=countMap1(:,i)/sum(countMap1(:,i));
end
for i=1:10
    plot(yBin1,countMap1(i,:),'Color',myColor(i,:,:));hold on
end
ends