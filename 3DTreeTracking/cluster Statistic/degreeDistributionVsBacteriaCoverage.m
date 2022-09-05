function [coverageTime,slopeTime]=degreeDistributionVsBacteriaCoverage(clusterTree,step,dirFile,bacteriaFrameInfo,matricType)
beginOne=step+1-mod(200,step)+200;
timeInterval=[beginOne:step:size(clusterTree,2),size(clusterTree,2)];
slopeTime=zeros(1,numel(timeInterval));
coverageTime=zeros(1,numel(timeInterval));
dirBioTreeResult=strcat(dirFile,'\bioTreeResult\t1\bioTree1.mat');
load(dirBioTreeResult)
imageProcessingInfo=bioTree{1}.imageProcessingInfo;
if strcmp(matricType,'w')
    for i=1:numel(timeInterval)
        iTime=timeInterval(i);
        degreeDistribution=getDistDistributionforWeighted(clusterTree{iTime}.weightedMatrix);
        if max(degreeDistribution(1,:))<=500
            slopeTime(i)=NaN;
        else
            slopeTime(i)=nicePolyFit(degreeDistribution(1,:),degreeDistribution(2,:),1);
        end
        coverageTime(i)=getCoverage(iTime,dirFile,imageProcessingInfo,bacteriaFrameInfo);
    end
else
    for i=1:numel(timeInterval)
        iTime=timeInterval(i);
        degreeDistribution=getDistDistribution(clusterTree{iTime}.distMatrix);
        if max(degreeDistribution(1,:))<=50
            slopeTime(i)=NaN;
        else
            slopeTime(i)=nicePolyFit(degreeDistribution(1,:),degreeDistribution(2,:),1);
        end
        coverageTime(i)=getCoverage(iTime,dirFile,imageProcessingInfo,bacteriaFrameInfo);
    end
end
end
function regionHist=getDistDistribution(distMatrix)
for i=1:size(distMatrix,2)
    iLine=distMatrix(i,:);
    degreeNum(i)=numel(iLine(iLine==1));
end
maxDegree=max(degreeNum);
regionHist(1,:)=[1:9,logspace(1,log10(maxDegree),30)];
for i=1:numel(regionHist(1,:))
    regionHist(2,i)=numel(degreeNum(degreeNum>=regionHist(1,i)));
end
end
function regionHist=getDistDistributionforWeighted(distMatrix)
for i=1:size(distMatrix,2)
    iLine=distMatrix(i,:);
    degreeNum(i)=sum(iLine(iLine~=0));
end
maxDegree=max(degreeNum);
regionHist(1,:)=[1:9,logspace(1,log10(maxDegree),60)];
for i=1:numel(regionHist(1,:))
    regionHist(2,i)=numel(degreeNum(degreeNum>=regionHist(1,i)));
end
end
function slope=nicePolyFit(xdata,ydata,beginD)
for i=10:numel(xdata)
    [~,b]=polyfit(log(xdata(beginD:i)),log(ydata(beginD:i)),1);
    if b.normr>0.05
        p=polyfit(log(xdata(beginD:i-1)),log(ydata(beginD:i-1)),1);
        slope=p(1)-1;
        hold on;plot(xdata,ydata)
        plot(xdata(beginD:i-1),exp(polyval(p,log(xdata(beginD:i-1)))),'r');
        break
    end
end
end
function coverageTime=getCoverage(iTime,dirFile,imageProcessingInfo,bacteriaFrameInfo)
% basicStackSize=189;
% iStack=fix((iTime-1)/basicStackSize)+1;
% nFrame=iTime-189*(iStack-1);
% tiffToMatFile=strcat(dirFile,'\tiff2matlab');
% nameList=dir(tiffToMatFile);
% properFile=strcat(tiffToMatFile,'\',nameList(iStack+2).name);
% load(properFile);
% imageStack=imageStack(:,:,nFrame);
% imageStack=myImageProcessingProper(imageStack,imageProcessingInfo);
% imageSize=size(imageStack);
% coverageTime=numel(imageStack(imageStack==1))/(imageSize(1)*imageSize(2));
cropInfo=imageProcessingInfo.cropInfo;
imageSize=numel(cropInfo(cropInfo==1));
bacteriaNum=size(bacteriaFrameInfo{iTime}.bacteriaInfo,1);
coverageTime=bacteriaNum/imageSize;
end
function afterProcessingImages=myImageProcessingProper(beforeProcessingImages,imageProcessingInfo) %this function can segreate orignal images as you want and returen the mask images
imageType='uint8'; %here you can chenge your image type
gaussianFilter=fspecial('gaussian',[5, 5],20); %here create Gaussian Blur Filter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
afterProcessingImages=zeros(size(beforeProcessingImages),imageType);
grayThresh=imageProcessingInfo.grayThresh;
areaThreshold=imageProcessingInfo.areaThreshold;
maxIntensity=imageProcessingInfo.maxIntensity;
for iframe=1:size(beforeProcessingImages,3)
    afterProcessingImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=im2bw(afterProcessingImages(:,:,iframe),grayThresh/255);%1.2 for the image with few particles
    afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MaxIntensity','PixelIdxList');
    for iCC=1:size(cc,1)
        if cc(iCC).MaxIntensity<=maxIntensity
            image=afterProcessingImages(:,:,iframe);
            image(cc(iCC).PixelIdxList)=0;
            afterProcessingImages(:,:,iframe)=image;
        end
    end
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'dilate'); % dilate process
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
end
afterProcessingImages=logical(afterProcessingImages);
end