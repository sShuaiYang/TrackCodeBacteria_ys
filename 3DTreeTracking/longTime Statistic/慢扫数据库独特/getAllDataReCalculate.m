function getAllDataReCalculate()
dirFile=uigetdir();
cd(dirFile)
nameList=dir(dirFile);
for i=3:numel(nameList)
    clc
    disp(nameList(i).name);
    name=nameList(i).name;
    if numel(name)>=3 && strcmp(name(1:3),'201')
        dirFinalTree=strcat(dirFile,'\',name,'\finalTree');
        cd(dirFinalTree)
        measureResult=dir(dirFinalTree);
        measureResultNum=numel(measureResult)-2;
        for j=1:measureResultNum
            load(strcat(dirFinalTree,'\',strcat('_',num2str(j),'_')));
            if j==1
                bioTree=bioTreeStack;
            else
                bioTree(end+1:end+size(bioTreeStack,2))=bioTreeStack;
            end
        end
        [bioTree,branchList,~,~]=myBiograph_new2(bioTree);
        bacteriaFrameInfo=getEachBacteriaInFrame(bioTree);
        allResult=dataAnalysisForLongTimeBacteria(bioTree);
        save(strcat(dirFile,'\',name,'\allResult'),'allResult');
        saveas(gcf,strcat(dirFile,'\',name,'\figure1.fig'));
        saveas(gcf,strcat(dirFile,'\',name,'\figure1.tif'));
        close all

        resultMap=traceStatisticMap(bioTree,bacteriaFrameInfo,20);
        imshow(resultMap)
        colormap(jet);
        saveas(gcf,strcat(dirFile,'\',name,'\traceInfo.tif'));
        close all
        
        load(strcat(dirFile,'\',name,'\allResult.mat'))
        [velocityResult,lengthTime]=lengthAndVelocityVsTime(bioTree,bacteriaFrameInfo);
        allResult.velocityResult=velocityResult;
        allResult.lengthTime=lengthTime;
        lengthAndVelocity(allResult.velocityResult,allResult.lengthTime);
        save(strcat(dirFile,'\',name,'\allResult'),'allResult');
        saveas(gcf,strcat(dirFile,'\',name,'\length&velocity.fig'))
        saveas(gcf,strcat(dirFile,'\',name,'\lengthAndvelocity.tif'),'tiffn')
%         allResult=dataAnalysisForLongTimeBacteria(bioTree);
%         save(strcat(dirFile,'\',name,'\allResult'),'allResult');
%         saveas(gcf,strcat(dirFile,'\',name,'\figure1.fig'))
        close all
        clear bioTree
        clear bioTreeStack
    end
end
end
function lengthAndVelocity(velocityResult,lengthTime)
load('\\192.168.1.163\d\ÂýÉ¨Êý¾Ý¿â\2014-02-11 F1 longTime(oilwater) Jzy\allResult.mat');
figure
hold on
subplot(1,2,2),plot(velocityResult(:,1),velocityResult(:,2),'r');
hold on;plot(allResult.velocityResult(:,1),allResult.velocityResult(:,2),'b');
xlabel('Time(min)');ylabel('velocity(um/s)');title('ave velocity vs Time');
subplot(1,2,1),plot(lengthTime(:,1),lengthTime(:,2),'r');
hold on;plot(allResult.lengthTime(:,1),allResult.lengthTime(:,2),'b');
xlabel('Time(min)');ylabel('length(um)');title('ave length vs Time');
end