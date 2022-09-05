function motilityDemoMaker(allData,startFrame,stepFrame,endFrame,dirFile)
% dirImages=uigetdir();
% dirSave=uigetdir();
dirImages=strcat(dirFile,'\tiff2matlab');
dirSave=strcat(dirFile,'\demo');
mkdir(dirSave)
figure;
h=gca;
cd(dirImages);
nameList=dir(dirImages);
vars = whos('-file', nameList(3).name);
stackSize=vars.size(3);
stackNum=fix(startFrame/(stackSize));
iImage=1;
if stackNum~=startFrame/(stackSize)
    fileName=nameList(stackNum+3).name;
else
    fileName=nameList(stackNum+2).name;
end
imageStack=loadImageStack(fileName);
for iframe=startFrame:stepFrame:endFrame;
    disp(iframe);
    stackNumNext=fix(iframe/(stackSize));
    if stackNumNext==iframe/(stackSize)
        stackNumNext=stackNumNext-1;
    end
    if stackNumNext~=stackNum
        clear imageStack;
        fileName=nameList(stackNumNext+3).name;
        imageStack=loadImageStack(fileName);
        stackNum=stackNumNext;
    end
    hold off;
    imshow(imageStack(:,:,iframe-stackNum*stackSize)); hold on; 
    plotDataPointandSave(dirSave,h,allData,iframe)
    iImage=iImage+1;
end
rmpath(dirImages);
end
function plotDataPointandSave(dirSave,h,allData,iframe)
    time_n=strcat(num2str(iframe),'frame');
    text(80,10,time_n,'color','g','FontSize',14);
    for iData=1:size(allData,2); % loop over all beads
%         label=getLabel(allData{iData});
        label=num2str(iData);
        [denoiseDataInFrame,velocityDataInFrame]=getDatainFrame(allData,iData,iframe);
        if ~isempty(denoiseDataInFrame)
            [p1VelocityVector,p2VelocityVector]=vectorPosition(denoiseDataInFrame,velocityDataInFrame);
            plot(denoiseDataInFrame(2),denoiseDataInFrame(3), ...
                'Marker','.','MarkerSize',3,'Color','g','LineStyle','none');hold on;
            plot(denoiseDataInFrame(4),denoiseDataInFrame(5), ...
                'Marker','.','MarkerSize',3,'Color','g','LineStyle','none');
            line([denoiseDataInFrame(2),denoiseDataInFrame(4)],[denoiseDataInFrame(3),denoiseDataInFrame(5)],'Marker','.','LineStyle','-','LineWidth',0.5)
            vectarrow(denoiseDataInFrame(2:3),p1VelocityVector);
            vectarrow(denoiseDataInFrame(4:5),p2VelocityVector);
            text(denoiseDataInFrame(2),denoiseDataInFrame(3),'1','color','g','FontSize',8);
            text(denoiseDataInFrame(4),denoiseDataInFrame(5),'2','color','g','FontSize',8);
            text(denoiseDataInFrame(6),denoiseDataInFrame(7),label,'color',[1,1,0],'FontSize',10);
        end
    end
    savefile2=strcat(dirSave,'\',num2str(iframe));
    saveas(h,savefile2,'tif');
end
function   [denoiseDataInFrame,velocityDataInFrame]=getDatainFrame(allData,iData,iframe)
denoiseData=allData{iData}.denosieData;
velocityData=allData{iData}.velocityData;
if iframe<denoiseData(1,1) || iframe>denoiseData(end,1)
    denoiseDataInFrame=[];
    velocityDataInFrame=[];
    return;
else
    frameIndex=iframe-denoiseData(1,1)+1;
    denoiseDataInFrame=denoiseData(frameIndex,:);
    velocityDataInFrame=velocityData(frameIndex,:);
end
end
function [p1VelocityVector,p2VelocityVector]=vectorPosition(denoiseDataInFrame,velocityDataInFrame)
p1p2Vector=denoiseDataInFrame(2:3)-denoiseDataInFrame(4:5);
p1p2Length=sqrt(sum(p1p2Vector.^2));
p1p2UnitVector=p1p2Vector./p1p2Length;
theat1=-velocityDataInFrame(5)*pi/180;
theat2=-velocityDataInFrame(10)*pi/180;
spinX1=p1p2UnitVector(1)*cos(theat1)-p1p2UnitVector(2)*sin(theat1);
spinY1=(p1p2UnitVector(1)*sin(theat1)+p1p2UnitVector(2)*cos(theat1));
spinX2=p1p2UnitVector(1)*cos(theat2)-p1p2UnitVector(2)*sin(theat2);
spinY2=(p1p2UnitVector(1)*sin(theat2)+p1p2UnitVector(2)*cos(theat2));
v1Amp=(log(velocityDataInFrame(4))+6)*3;
v2Amp=(log(velocityDataInFrame(9))+6)*3;
p1VelocityVector=[denoiseDataInFrame(2)+spinX1*v1Amp,denoiseDataInFrame(3)+spinY1*v1Amp];
p2VelocityVector=[denoiseDataInFrame(4)+spinX2*v2Amp,denoiseDataInFrame(5)+spinY2*v2Amp];
end
function vectarrow(p0,p1)
if max(size(p1))==2
    x0 = p0(1);
    y0 = p0(2);
    x1 = p1(1);
    y1 = p1(2);
    plot([x0;x1],[y0;y1],'Color',[1 0 0],'LineStyle','-','LineWidth',1);   % Draw a line between p0 and p1
    p = p1-p0;
    alpha = 0.5;  % Size of arrow head relative to the length of the vector
    beta = 0.5;  % Width of the base of the arrow head relative to the length
    hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
    hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
    hold on
    plot(hu(:),hv(:),'Color',[1 0 0],'LineStyle','-','LineWidth',1)  % Plot arrow head
end
end
function imageStack=loadImageStack(fileName)
tempImage=load(fileName);
imageStack=tempImage.imageStack;
end
function label=getLabel(dataList)
if ~isempty(dataList.branchIndex)
    branchIndex=dataList.branchIndex;
    if dataList.isNode==false
        rootInfo=dataList.rootInfo;
        label=strcat('B',num2str(branchIndex),'R',num2str(rootInfo(1)),'/',num2str(rootInfo(2))) ;
    elseif dataList.isNode==true
        nodeInfo=dataList.nodeInfo;
        label=strcat('B',num2str(branchIndex),'N',num2str(nodeInfo(1)),'/',num2str(nodeInfo(2)),'/',num2str(nodeInfo(3))) ;
    end
elseif isempty(dataList.branchIndex)
    if dataList.isNode==false
        rootInfo=dataList.rootInfo;
        label=strcat('R',num2str(rootInfo(1)),'/',num2str(rootInfo(2))) ;
    end
end
end
