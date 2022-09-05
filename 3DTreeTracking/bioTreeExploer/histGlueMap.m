function histResult=histGlueMap(timeGlueMap)
for itime=1:size(timeGlueMap,3)
    I=timeGlueMap(:,:,itime);
    histTemp=zeros(1,max(max(I)));
    for iNum=1:max(max(I))
        histTemp(iNum)=size(I(I==iNum),1);
    end
    histResult(itime).data=histTemp;
end
end