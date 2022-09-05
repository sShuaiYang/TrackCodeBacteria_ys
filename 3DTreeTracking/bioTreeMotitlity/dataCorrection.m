function allData=dataCorrection(allDataOri,correctionList,referenceNum,ZDriftRef,dirSave,rockFrame)
% dirSave=uigetdir();
dirSave1=strcat(dirSave,'\','allDataOri');
save(dirSave1,'allDataOri');
i=0;
%here we choose the data we need,because some picture may have influence in
%our result. Then we inverse the head and tail of some moving
%bacterias. '1' will be the head. Another change is that we have the
%variable allDataOri type 'cell' to a structure type allData. It has three
%parts.DataInfo for information of different bacterias,correctionList for
%the information we used that cha nge the allDataOri,dataList tells us the
%new number we give for each bacteria.
for iData=1:size(allDataOri,2)
    dataTemp=allDataOri{iData};
    frameCutAndInverseInfo=correctionList(correctionList(:,1)==iData,2:end);
    if ~isempty(frameCutAndInverseInfo)
        for small_iData=1:ceil(size(frameCutAndInverseInfo,2)/3)
            frameCut=frameCutAndInverseInfo(1,3*small_iData-2:3*small_iData);
            if ~isempty(frameCut(frameCut>0))
                i=i+1;
                dataList{i}=[iData,small_iData];
                allData.dataInfo{i}=dataTemp;
                allData.dataInfo{i}.originalData=dataTemp.originalData(dataTemp.originalData(:,1)>=frameCut(1) & dataTemp.originalData(:,1)<=frameCut(2),:);
                allData.dataInfo{i}.denosieData=dataTemp.denosieData(dataTemp.originalData(:,1)>=frameCut(1) & dataTemp.originalData(:,1)<=frameCut(2),:);
                allData.dataInfo{i}.velocityData=dataTemp.velocityData(dataTemp.originalData(:,1)>=frameCut(1) & dataTemp.originalData(:,1)<=frameCut(2),:);  %select data we need
                if frameCut(3)==1
                    allData=inverseData(allData,i);     %inverse the head and tail
                end
                allData.dataInfo{i}.denosieData(:,9)=allData.dataInfo{i}.denosieData(:,8);   %denosieData(:,9)is the major axis of a bacteria,the same as original Data
                allData.dataInfo{i}.denosieData(:,8)=findAngle([allData.dataInfo{i}.denosieData(:,2)-allData.dataInfo{i}.denosieData(:,4),-(allData.dataInfo{i}.denosieData(:,3)-allData.dataInfo{i}.denosieData(:,5))]);
                %denoisedata(:,8) is the angle between the bacteria major axis and X axis
            end
        end
    end
end
allData.rockFrame=rockFrame;
allData.correctionList=correctionList;
allData.dataList=dataList;
allData.XYdriftRef=referenceNum;
allData.ZDrift=ZDriftRef;
%here we correct the Centroid position
%we must choose a bacteria that is nearly motionless,to have the Centroid
%postion correction.If we can't find this becteria,the referenceNum=[];for
%one bacteria,if the picture isn't continuous,we can choose different parts
%of this bacteria,like[1.1,1.2,2.1]for one bacteria.

%Attention:we'd best choose only one bacteria in this function,our function
%          is just believable for only one bacteria
if ~isempty(referenceNum)
    allData=centroidCorrection(allData,referenceNum);
else
    disp('No XY Drift Correction(centroid Correcition)')
end

%here we correct the length and Oritation
allData=allDataCheck(allData,ZDriftRef);
saceFile2=strcat(dirSave,'\allData');
save(saceFile2,'allData');
end

%this function is used for inverse the sequence of a bacteria,show which is
%its 'head',because we will caculate the velocity information after we correct sita, lens and centroid position, as a result, we ignore the velocity here.
function allData=inverseData(allData,i)
a=allData.dataInfo{i}.originalData(:,2:3);
allData.dataInfo{i}.originalData(:,2:3)=allData.dataInfo{i}.originalData(:,4:5);
allData.dataInfo{i}.originalData(:,4:5)=a;

a=allData.dataInfo{i}.denosieData(:,2:3);
allData.dataInfo{i}.denosieData(:,2:3)=allData.dataInfo{i}.denosieData(:,4:5);
allData.dataInfo{i}.denosieData(:,4:5)=a;
end

function angle=findAngle(vector) %return the angle ([-180,180]) betwwen the vector1 and axis X

angle=atan(vector(:,2)./abs(vector(:,1))).*(180/pi);
angle(vector(:,1)<0&vector(:,2)>=0)= 180-angle(vector(:,1)<0&vector(:,2)>=0);
angle(vector(:,1)<0&vector(:,2)<0)= -180-angle(vector(:,1)<0&vector(:,2)<0);

angle(angle>180)=angle(angle>180)-360;
angle(angle<-180)=angle(angle<-180)+360;
end

%At first, we combine all the time series information of a bacteria; then
%we caculate the x-y drift;at last every bacteria will minus this drift.
function allData=centroidCorrection(allData,referenceNum)
for i=1:numel(allData.dataList)
    [p,~]=find(referenceNum==allData.dataList{i}(1)+0.1*allData.dataList{i}(2));
    if ~isempty(p)
        referenceBecteriaCentroidX(allData.dataInfo{i}.denosieData(:,1),1)=allData.dataInfo{i}.denosieData(:,6);
        referenceBecteriaCentroidY(allData.dataInfo{i}.denosieData(:,1),1)=allData.dataInfo{i}.denosieData(:,7);
    end
end
referenceBecteriaCentroidX=fillIn(referenceBecteriaCentroidX);
referenceBecteriaCentroidY=fillIn(referenceBecteriaCentroidY);
diffXData=getDiffData(referenceBecteriaCentroidX);
diffYData=getDiffData(referenceBecteriaCentroidY);
timeData=1:numel(diffXData);
for iData=1:size(allData.dataInfo,2)
    timeInfo=allData.dataInfo{iData}.denosieData(:,1);
    if timeInfo(end)>timeData(end)
        diffXData(timeData(end)+1:timeInfo(end),1)=diffXData(timeData(end));
        diffYData(timeData(end)+1:timeInfo(end),1)=diffYData(timeData(end));
        timeData=1:numel(diffXData);
    end
    allData.dataInfo{iData}.denosieData(:,6)=allData.dataInfo{iData}.denosieData(:,6)-diffXData(timeInfo(timeInfo>=timeData(1)&timeInfo<=timeData(end)));
    allData.dataInfo{iData}.denosieData(:,6)=sectionDenoise(allData.dataInfo{iData}.denosieData(:,6),0.2,7);
    allData.dataInfo{iData}.denosieData(:,7)=allData.dataInfo{iData}.denosieData(:,7)-diffYData(timeInfo(timeInfo>=timeData(1)&timeInfo<=timeData(end)));
    allData.dataInfo{iData}.denosieData(:,7)=sectionDenoise(allData.dataInfo{iData}.denosieData(:,7),0.2,7);
end
end

%here we have the information of the motionless bacteria and regard its
%track as the x-y drift of the platform.So we smooth it and caculate how much it
%changes.
function diffData=getDiffData(oriData)
[~,sorh,keepapp] = ddencmp('den','wv',oriData);
thr=1/0.01;
NewData=wdencmp('gbl',oriData,'db7',7,thr,sorh,keepapp);
diffData=NewData-oriData(1);
end

function data=fillIn(data)
if data(1)==0
    [p,~]=find(data~=0);
    data(1)=data(p(1));
end
for i=1:numel(data)
    if data(i)==0
        data(i)=data(i-1);
    end
end
end

%wdencmp we used will decrease the noise,but also remove some useful
%information we need. So we use this method for different sections.For
%those data whose change below the threshold, we use the wdencmp function,
%for others,we keep the original data.
function newSignal=sectionDenoise(originalSignal,threshold,strength)
diffSignal=diff(originalSignal);
diffSignal=[0;diffSignal];
position=find(abs(diffSignal)>=threshold);
position=[1;position;numel(diffSignal)];
for i=1:numel(position)-2
    Data=originalSignal(position(i):position(i+1)-1,1);
    if numel(Data)>=3
        [~,sorh,keepapp] = ddencmp('den','wv',Data);
        thr=1/0.01;
        newSignal(position(i):position(i+1)-1,1)=wdencmp('gbl',Data,'db7',strength,thr,sorh,keepapp);
    else
        newSignal(position(i):position(i+1)-1,1)=Data;
    end
end
Data=originalSignal(position(end-1):position(end),1);
[~,sorh,keepapp] = ddencmp('den','wv',Data);
thr=1/0.01;
newSignal(position(end-1):position(end),1)=wdencmp('gbl',Data,'db7',strength,thr,sorh,keepapp);
end