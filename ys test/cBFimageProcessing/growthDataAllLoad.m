function  growthDataAll=growthDataAllLoad(dirFile)
% growthDataAllLoad
% dirFile='E:\2019-12-23 PAO1_IP32_100x ys\growthDataAll';
nameList=dir(dirFile);
for i=1:numel(nameList)-2
    load( [dirFile,'\',nameList(i+2).name]);
    if i==1
        growthDataAll=growthDataSub;
    else
        growthDataAll(end+1:end+numel(growthDataSub))=growthDataSub;
    end       
end

end