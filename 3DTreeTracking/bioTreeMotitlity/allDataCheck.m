function allData=allDataCheck(allData,ZDriftRef)
%this function has two parts:bacteria length correction and sita correction

%here we have the bacteria length correction.The change of bacteria length include three parts:growing up of a bacteria,Z drift and whether it will moving via standing up.
%At first,we choose some bacteria that will not move at Z positon(moving
%without standing up). Then we make a linearization for each bacteria we
%choose and caculate the Z drift.At last all the bacteria length will minus
%this value.
if ~isempty(ZDriftRef)
for iData=1:numel(allData.dataList)
    cell=allData.dataList{iData};
    [u,~]=find(ZDriftRef==cell(1)+0.1*cell(2));
    if ~isempty(u)
        time=allData.dataInfo{iData}.denosieData(:,1);
        bacteriaLen(time,cell(1))=sqrt((allData.dataInfo{iData}.denosieData(:,2)-allData.dataInfo{iData}.denosieData(:,4)).^2+(allData.dataInfo{iData}.denosieData(:,3)-allData.dataInfo{iData}.denosieData(:,5)).^2);
        maxY=max(bacteriaLen(time,cell(1)));
        changeY=maxY-bacteriaLen(time,cell(1));
        NewX=time(changeY<=4e-4*numel(time));
        t=bacteriaLen(time,cell(1));
        NewY=t(changeY<=4e-4*numel(time));
        p=polyfit(NewX,NewY,1);
        yData=polyval(p,time);
        driftData(time,cell(1))=bacteriaLen(time,cell(1))-yData;
    end
end
for i=1:size(driftData,1)
    rowData=driftData(i,:);
    driftDataAll(i,1)=mean(rowData(rowData~=0));
end
driftDataAll(isnan(driftDataAll))=0;
driftDataAll=fillIn(driftDataAll);
[~,sorh,keepapp] = ddencmp('den','wv',driftDataAll);
thr=1/0.01;
driftDataAll=wdencmp('gbl',driftDataAll,'db7',6,thr,sorh,keepapp);
clear bacteriaLen
else 
    driftDataAll=[];
    disp('No Z Drift Correction(length Correcition)')
end
for iData=1:numel(allData.dataList)
    allData.dataInfo{iData}.DenoisedLenInfo=sqrt((allData.dataInfo{iData}.denosieData(:,2)-allData.dataInfo{iData}.denosieData(:,4)).^2+(allData.dataInfo{iData}.denosieData(:,3)-allData.dataInfo{iData}.denosieData(:,5)).^2);
    if size(driftDataAll,1)<max(allData.dataInfo{iData}.denosieData(:,1))
        driftDataAll(end+1:max(allData.dataInfo{iData}.denosieData(:,1)),1)=0;
    end
    allData.dataInfo{1,iData}.beforeWdencmpLen=allData.dataInfo{iData}.DenoisedLenInfo-driftDataAll(allData.dataInfo{iData}.denosieData(:,1));
    allData.dataInfo{1,iData}.afterWdencmpLen=sectionDenoise(allData.dataInfo{1,iData}.beforeWdencmpLen,0.5,7);
end


%here we have the Oritation correction.So here we have already the
%sita,length and centroid positon.We can gain the head and the tail of the
%bacteria easily.Then we can caculate all the velocity information we need.
for iData=1:numel(allData.dataInfo)
    allData.dataInfo{iData}.denosieData(:,8)=getOritationbywdencmp(allData.dataInfo{iData}.denosieData(:,8));
    [allData.dataInfo{iData}.denosieData(:,2:3),allData.dataInfo{iData}.denosieData(:,4:5)]=twoPointPosition(allData.dataInfo{iData}.denosieData(:,6:7),allData.dataInfo{iData}.denosieData(:,8),allData.dataInfo{iData}.afterWdencmpLen,1);
    [~,allData.dataInfo{iData}.velocityData(1:end-1,2:16)]=getVelocity(allData.dataInfo{iData}.denosieData(:,[2:7,9]));
    allData.dataInfo{iData}.velocityData(end,:)=[];
    allData.dataInfo{iData}.denosieData(end,:)=[];
end
end

%this function is used to gain the two polars points of a
%bacteria.Eccentricity here we use 1,because the length of bacteria here is
%actually 2*c.
function [p1Position,p2Position]=twoPointPosition(centroid,orientation,majorAxisLength,eccentricity)
O_sign=sign(orientation);
O_cos=cos(abs(pi*orientation/180));
O_sin=sin(abs(pi*orientation/180));
Axis_length=0.5*majorAxisLength;
x=centroid(:,1);
y=centroid(:,2);
p1_x=x+Axis_length.*O_cos*eccentricity;
p1_y=y-O_sign.*Axis_length.*O_sin*eccentricity;
p2_x=x-Axis_length.*O_cos.*eccentricity;
p2_y=y+O_sign.*Axis_length.*O_sin*eccentricity;
p1Position=[p1_x,p1_y];
p2Position=[p2_x,p2_y];
end

%dataList has 7 rows. row[1:2](1x 1y), row[3:4](2x 2y), row[5:6](Cx Cy),
%row[7] is the length of the bacteria.
function [dataListDenoise,dataVelocity] =getVelocity(dataList)
dataListDenoise=dataList;
if size(dataListDenoise,1)>1
    dataVelocityTemp=zeros(size(dataList,1)-1,size(dataList,2));
    for iData=1:size(dataList,2)
        dataVelocityTemp(:,iData)=diff(dataListDenoise(:,iData));
    end
    p1VelocityAmp=sqrt(dataVelocityTemp(:,1).^2+dataVelocityTemp(:,2).^2);
    p2VelocityAmp=sqrt(dataVelocityTemp(:,3).^2+dataVelocityTemp(:,4).^2);
    centroidVelocityAmp=sqrt(dataVelocityTemp(:,5).^2+dataVelocityTemp(:,6).^2);
    angleV1Pre=angleBetweenTwoVector(dataVelocityTemp(:,1:2),[dataListDenoise(1:end-1,1)-dataListDenoise(1:end-1,3),dataListDenoise(1:end-1,2)-dataListDenoise(1:end-1,4)]);
    angleV2Pre=angleBetweenTwoVector(dataVelocityTemp(:,3:4),[dataListDenoise(1:end-1,1)-dataListDenoise(1:end-1,3),dataListDenoise(1:end-1,2)-dataListDenoise(1:end-1,4)]);
    angleV1next=angleBetweenTwoVector(dataVelocityTemp(:,1:2),[dataListDenoise(2:end,1)-dataListDenoise(2:end,3),dataListDenoise(2:end,2)-dataListDenoise(2:end,4)]);
    angleV2next=angleBetweenTwoVector(dataVelocityTemp(:,3:4),[dataListDenoise(2:end,1)-dataListDenoise(2:end,3),dataListDenoise(2:end,2)-dataListDenoise(2:end,4)]);
    angleVelocity=angleBetweenTwoVector([dataListDenoise(2:end,1)-dataListDenoise(2:end,3),dataListDenoise(2:end,2)-dataListDenoise(2:end,4)],[dataListDenoise(1:end-1,1)-dataListDenoise(1:end-1,3),dataListDenoise(1:end-1,2)-dataListDenoise(1:end-1,4)]);
    dataVelocity(:,1:2)=dataVelocityTemp(:,1:2);
    dataVelocity(:,3)=p1VelocityAmp;
    dataVelocity(:,4)= angleV1Pre;
    dataVelocity(:,5)= angleV1next;
    dataVelocity(:,6:7)=dataVelocityTemp(:,3:4);
    dataVelocity(:,8)=p2VelocityAmp;
    dataVelocity(:,9)= angleV2Pre;
    dataVelocity(:,10)= angleV2next;
    dataVelocity(:,11:12)=dataVelocityTemp(:,5:6);
    dataVelocity(:,13)=centroidVelocityAmp;
    dataVelocity(:,14)=angleVelocity;
    dataVelocity(:,15)=dataVelocityTemp(:,7);
else
    dataVelocity=[];
end
end

function angle=angleBetweenTwoVector(vector1,vector2) %return the angle ([-180,180]) betwwen the vector1 and vector2

angle1=atan(vector1(:,2)./abs(vector1(:,1))).*(180/pi);
angle1(vector1(:,1)<0&vector1(:,2)>=0)= 180-angle1(vector1(:,1)<0&vector1(:,2)>=0);
angle1(vector1(:,1)<0&vector1(:,2)<0)= -180-angle1(vector1(:,1)<0&vector1(:,2)<0);

angle2=atan(vector2(:,2)./abs(vector2(:,1))).*(180/pi);
angle2(vector2(:,1)<0&vector2(:,2)>=0)= 180-angle2(vector2(:,1)<0&vector2(:,2)>=0);
angle2(vector2(:,1)<0&vector2(:,2)<0)= -180-angle2(vector2(:,1)<0&vector2(:,2)<0);

angle=angle2- angle1;
angle(angle>180)=angle(angle>180)-360;
angle(angle<-180)=angle(angle<-180)+360;
end

%correction for sita,we just need to smooth it to remove the noise.Here we
%use angleDisplacement instead of Oritation,because we consider that
%bacteria will not rotate above 180 degree.
function NewOritation=getOritationbywdencmp(angleToTime)
angleVelocity=diff(angleToTime);
angleVelocity(angleVelocity>180)=angleVelocity(angleVelocity>180)-360;
angleVelocity(angleVelocity<-180)=angleVelocity(angleVelocity<-180)+360;
angleDisplacement=zeros(numel(angleToTime),1);
for i=1:numel(angleToTime)
    angleDisplacement(i)=sum(angleVelocity(1:i-1))+angleToTime(1);
end
NewOritation=sectionDenoise(angleDisplacement,1,7);
end

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