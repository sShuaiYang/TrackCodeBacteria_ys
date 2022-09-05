function [glueMapC,allRg,rgList]=glueMapExplorer_new(glueMap,glueMapSearching)
rgMassList=[];
rgList=[];
xSize=size(glueMap,1);
ySize=size(glueMap,2);
frameNum=size(glueMap,3);
glueMapC=false(xSize,ySize,frameNum);
allRg=zeros(1,frameNum);
allRgs=zeros(1,frameNum);
broderWithe=5; % here create broder frame
broderFrame=true(xSize,ySize);broderFrame(:,1:broderWithe)=0;broderFrame(:,ySize-broderWithe:ySize)=0;broderFrame(1:broderWithe,:)=0;broderFrame(xSize-broderWithe:xSize,:)=0;
parfor iframe=1:frameNum
%     dispFrame(iframe)
    glueMapC(:,:,iframe)=~glueMap(:,:,iframe)&broderFrame;
    CC=bwconncomp(glueMapC(:,:,iframe));
    stats = regionprops(CC, 'Centroid','Area','PixelList');
    if ~isempty(stats)
        [Rg,~]=averageRg(stats);
        allRg(1,iframe)=Rg;
%         rgMassList=[rgMassList,rgListTemp(:,2)'];
%         rgList=[rgList,rgListTemp(:,1)'];
    end
end
parfor iframe=1:frameNum
%     dispFrame(iframe)
    CC=bwconncomp(glueMapSearching(:,:,iframe));
    stats = regionprops(CC, 'Centroid','FilledArea','PixelList');
    if ~isempty(stats)
        [Rgs,rgListTemp]=averageRgs(stats);
        allRgs(1,iframe)=Rgs;
        rgMassList=[rgMassList,rgListTemp(:,2)'];
        rgList=[rgList,rgListTemp(:,1)'];
    end
end
rgList=[rgList',rgMassList'];
allRg=[allRg',allRgs'];
end
function [Rg,rgList]=averageRg(stats)
rgList=zeros(size(stats,1),2);
for iObjective=1:size(stats,1)
    X0=stats(iObjective).Centroid(1);
    Y0=stats(iObjective).Centroid(2);
    rgList(iObjective,2)=stats(iObjective).Area;
    rgList(iObjective,1)=sqrt(mean((stats(iObjective).PixelList(:,1)-X0).^2+(stats(iObjective).PixelList(:,2)-Y0).^2));
end
Rg=sum(rgList(:,1).*rgList(:,2))/sum(rgList(:,2));
end
function [Rg,rgList]=averageRgs(stats)
rgList=zeros(size(stats,1),2);
for iObjective=1:size(stats,1)
    X0=stats(iObjective).Centroid(1);
    Y0=stats(iObjective).Centroid(2);
    rgList(iObjective,2)=stats(iObjective).FilledArea;
    rgList(iObjective,1)=sqrt(mean((stats(iObjective).PixelList(:,1)-X0).^2+(stats(iObjective).PixelList(:,2)-Y0).^2));
end
Rg=sum(rgList(:,1).*rgList(:,2))/sum(rgList(:,2));
end
function dispFrame(iConnect) 
lengthB=length(num2str(iConnect-1));
back='';
for i=1:lengthB
    back=strcat(back,'\b');
end
fprintf(strcat(back,num2str(iConnect)));
end