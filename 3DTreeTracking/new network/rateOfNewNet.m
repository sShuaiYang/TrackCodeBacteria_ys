function resultAll=rateOfNewNet()
dirFile=uigetdir();
nameList=dir(dirFile);
n=0;
for i=1:numel(nameList)-2
    nameI=nameList(i+2).name;
    dirFile1=strcat(dirFile,'\',nameI);
    cd(dirFile1);
    nameList1=dir(dirFile1);
    for iResult=1:numel(nameList1)-2
        load(strcat(dirFile1,'\',num2str(iResult),'.mat'))
        for iStack=1:size(result,2)
            n=n+1;
            disp(n)
            if iResult==1 && i==1 && iStack==1
                for iSmall=1:24
                    resultAll{iSmall}.frameNum=200*iSmall;
                    resultAll{iSmall}.finalMatrix=result{iStack}(iSmall).finalMatrix;
                end
                continue
            end
            for iSmall=1:24
                resultAll{iSmall}.finalMatrix=result{iStack}(iSmall).finalMatrix+resultAll{iSmall}.finalMatrix;
            end
        end
    end       
end
end