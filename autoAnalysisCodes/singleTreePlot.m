function singleTreePlot(allData)
% allData里面含两类数据，一个是从mask得到的bacInfo
% 第一列为time(min),第二列为majorLength,第三列为minorLength,第四列为CyOFP,第五列为GFP,第六列为mScalet,第七列为RFP,第八列为CyPet,第九列为Venus,第十列为Ametrine一一对应，没有数据的地方用NaN表示
% 第一列为time(min),第二列为majorLength,第三列为minorLength,第四列为CyOFP,第五列为GFP,第六列为mScalet,第七列为RFP，一一对应，没有数据的地方用NaN表示，第八列为Orientation,第九列为FilledArea,
% 第十列为growthRate,第十一列为fpProduceGFP，第十二列为treeSize,第十三列为generation,第十四列为CyPet,第十五列为Venus,第十六列为Ametrine，第十七列为fpProduceCyOFP,第十八列为fpProducemScalet
% 第十九列为fpProduceRFP，第二十列为fpProduceCyPet，第二十一列为fpProduceVenus,第二十二列为fpProduceAmetrine,第二十三列为redControl,第二十四列为blueControl,第二十五列为greenControl
clc
iNum=numel(allData.treeAll);
disp(['please input tree num (1 to ',num2str(iNum),')']);
treeNum=input('____:');
disp('please input treePlot type')
variableName=checkData(allData);
treeType=input('____:');
treeType=variableName{treeType};
try
    drawDataMap(allData.treeAll,treeNum,treeType);
catch err
    disp('do not have this fluo color')
end
end
function [xData,yData]=drawDataMap(allData,varargin)
linkMatrix=zeros(size(allData(varargin{1}).linkMatrix));
linkRow=allData(varargin{1}).linkRow;
linkColumn=allData(varargin{1}).linkColumn;
for i=1:numel(linkRow)
    switch varargin{2}
        case 'CyOFP'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).CyOFP(~isnan(allData(varargin{1}).linkInfo(i).CyOFP)));
        case 'GFP'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).GFP(~isnan(allData(varargin{1}).linkInfo(i).GFP)));
        case 'mScalet'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).mScalet(~isnan(allData(varargin{1}).linkInfo(i).mScalet)));
        case 'RFP'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).RFP(~isnan(allData(varargin{1}).linkInfo(i).RFP)));
        case 'CyPet'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).CyPet(~isnan(allData(varargin{1}).linkInfo(i).CyPet)));
        case 'Venus'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).Venus(~isnan(allData(varargin{1}).linkInfo(i).Venus)));
        case 'mAmetrine'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).mAmetrine(~isnan(allData(varargin{1}).linkInfo(i).mAmetrine)));
        case 'Age'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).linkAge);
        case 'simpleTree'
            linkMatrix(linkRow(i),linkColumn(i))=20;
        case 'time'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).linkTimer);
        case 'maxLen'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).measurment{ceil(allData(varargin{1}).linkInfo(i).measurment/2)}.MajorAxisLength);
        case 'minLen'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).measurment{ceil(allData(varargin{1}).linkInfo(i).measurment/2)}.MinorAxisLength);
        case 'growthRate'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).growthRate(~isnan(allData(varargin{1}).linkInfo(i).growthRate)));
        case 'fpProduceCyOFP'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).fpProduceCyOFP(~isnan(allData(varargin{1}).linkInfo(i).fpProduceCyOFP)));
        case 'fpProduceGFP'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).fpProduceGFP(~isnan(allData(varargin{1}).linkInfo(i).fpProduceGFP)));
        case 'fpProducemScalet'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).fpProducemScalet(~isnan(allData(varargin{1}).linkInfo(i).fpProducemScalet)));
        case 'fpProduceRFP'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).fpProduceRFP(~isnan(allData(varargin{1}).linkInfo(i).fpProduceRFP)));
        case 'fpProduceCyPet'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).fpProduceCyPet(~isnan(allData(varargin{1}).linkInfo(i).fpProduceCyPet)));
        case 'fpProduceVenus'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).fpProduceVenus(~isnan(allData(varargin{1}).linkInfo(i).fpProduceVenus)));
        case 'fpProducemAmetrine'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).fpProducemAmetrine(~isnan(allData(varargin{1}).linkInfo(i).fpProducemAmetrine)));
        case 'redControl'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).redControl);
        case 'blueControl'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).blueControl);
        case 'greenControl'
            linkMatrix(linkRow(i),linkColumn(i))=mean(allData(varargin{1}).linkInfo(i).greenControl);
    end
end
plotLinkTree(linkMatrix,allData(varargin{1}).timer,1:allData(varargin{1}).leafNum,allData(varargin{1}).cutNum);
hold on
ylabel(['Field',num2str(allData.fieldNum),':',varargin{2}])
end
function plotLinkTree(linkMatrix,centroidInfo,leafYCoo,cutNum)
% generally linkMatrix means the relationship for [N, L], N means all
% node, L means all leaf, centroidInfo here could be changed for the second
% column. The first column means the frame that the node begin\
% color=[1,0,1];
leafNum=numel(leafYCoo);
nodeNum=size(linkMatrix,1)-leafNum;
centroidInfo(nodeNum+1:end,2)=leafYCoo;
leafOrder=1:leafNum;
for i=1:nodeNum
    iLine=linkMatrix(nodeNum+1-i,nodeNum+1-i:end);
    centroid2=centroidInfo(:,2);
    centroid2=centroid2(nodeNum+1-i:end);
    leafPos=(centroid2(iLine~=0));
    centroidInfo(nodeNum+1-i,2)=mean(leafPos(~isnan(leafPos)));
end
% figure;
colorAll=colormap(jet(124));
colorAll=colorAll(25:124,:);
close all
maxLinkNum=max(max(linkMatrix(linkMatrix~=0)));
minLinkNum=min(min(linkMatrix(linkMatrix~=0)));
colorRange=input(['please tap the color range ',num2str(minLinkNum) ,' to ',num2str(maxLinkNum),':']);
minLinkNum=colorRange(1);
maxLinkNum=colorRange(2);

% new add map size
figure1=figure;
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.288229166666667 0.74854189336235],'FontSize',20);
% ylim(axes1,[0 118])
% xlim(axes1,[0 12.6])
view(axes1,[90 -90]);
box(axes1,'on');
hold all
% Create xlabel
xlabel('time(h)','VerticalAlignment','bottom','Rotation',90,...
    'HorizontalAlignment','center','FontSize',20);
scrsz=get(0,'ScreenSize');
set (gcf,'Position',[10 10 scrsz(3) scrsz(4)-80]);
set(gcf, 'PaperPositionMode', 'auto');

for i=1:size(linkMatrix,1)
    for j=i+1:size(linkMatrix,1)
        if linkMatrix(i,j)~=0
            if maxLinkNum==minLinkNum
                colorIndex=1;
            else
                colorIndex=round(99*(linkMatrix(i,j)-minLinkNum)/(maxLinkNum-minLinkNum))+1;
            end
            if colorIndex<=1
                colorIndex=1;
            end
            if colorIndex>=100;
                colorIndex=100;
            end
            linkTwoPointsAging(centroidInfo(i,:),centroidInfo(j,:),nodeNum,i,j,cutNum,colorAll(colorIndex,:));
        end
    end
end
end
function linkTwoPointsAging(c1,c2,nodeNum,i,j,cutNum,colorI)
c1(1)=c1(1)/60;
c2(1)=c2(1)/60;
cutNum=cutNum/60;
if c1(1)>c2(1)
    c=c2;
    c2=c1;
    c1=c;
end
hold on;
if c1(1)>cutNum
    return
end
% link line
% c3=[c1(1),c2(2)];
% line([c1(1),c3(1),c2(1)],[c1(2),c3(2),c2(2)],'Color',color,'LineWidth',1);

% link with ellipse
c3=[c2(1),c1(2)];
axisA=abs(c2(1)-c1(1));
axisB=abs(c2(2)-c1(2));
if c2(1)>cutNum
    x=c1(1):0.01:cutNum;
else
    x=c1(1):0.01:c2(1);
end
if c2(2)<c1(2)
    y=c3(2)-axisB*(1-((x-c3(1))/axisA).^2).^(1/2);
else
    y=c3(2)+axisB*(1-((x-c3(1))/axisA).^2).^(1/2);
end
line(x,y,'Color',colorI,'LineWidth',1)

maker='.';
size=16;
if i==1
    plot(c1(1),c1(2),'Marker','^','MarkerSize',12,'Color',[1,0,0]);
else
    if j<=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
    end
    if j>=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        if c2(1)<360
            plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[85,107,47]/255);
        end
    end
end
end
function variableName=checkData(allData)
% allData里面含两类数据，一个是从mask得到的bacInfo
% 第一列为time(min),第二列为majorLength,第三列为minorLength,第四列为CyOFP,第五列为GFP,第六列为mScalet,第七列为RFP,第八列为CyPet,第九列为Venus,第十列为Ametrine一一对应，没有数据的地方用NaN表示
% 第一列为time(min),第二列为majorLength,第三列为minorLength,第四列为CyOFP,第五列为GFP,第六列为mScalet,第七列为RFP，一一对应，没有数据的地方用NaN表示，第八列为Orientation,第九列为FilledArea,
% 第十列为growthRate,第十一列为fpProduceGFP，第十二列为treeSize,第十三列为generation,第十四列为CyPet,第十五列为Venus,第十六列为Ametrine，第十七列为fpProduceCyOFP,第十八列为fpProducemScalet
% 第十九列为fpProduceRFP，第二十列为fpProduceCyPet，第二十一列为fpProduceVenus,第二十二列为fpProduceAmetrine,第二十三列为redControl,第二十四列为blueControl,第二十五列为greenControl
variableName{1}='time';
variableName{2}='maxLen';
variableName{3}='minLen';
variableName{4}='CyOFP';
variableName{5}='GFP';
variableName{6}='mScalet';
variableName{7}='RFP';
variableName{8}='CyPet';
variableName{9}='Venus';
variableName{10}='mAmetrine';
variableName{11}='angle';
variableName{12}='area';
variableName{13}='growthRate';
variableName{14}='tag';
variableName{15}='tagValue';
variableName{16}='Age';
variableName{17}='generation';
variableName{18}='fpProduceCyOFP';
variableName{19}='fpProduceGFP';
variableName{20}='fpProducemScalet';
variableName{21}='fpProduceRFP';
variableName{22}='fpProduceCyPet';
variableName{23}='fpProduceVenus';
variableName{24}='fpProducemAmetrine';
variableName{25}='redControl';
variableName{26}='blueControl';
variableName{27}='greenControl';
variableName{28}='simpleTree';
maskImage=allData.maskResult;
trackingResult=allData.trackingResult;
chooseIndex=true(28,1);
for i=4:10
    if all(isnan(maskImage(:,i))) 
        chooseIndex(i)=0;
        chooseIndex(i+14)=0;
    end
end
for i=25:27
    if all((trackingResult(:,i-2)==0))
        chooseIndex(i)=0;
    end
end
variableName=variableName(chooseIndex);
variableNum=numel(variableName);
line1=[];
line2=[];
line3=[];
for i=1:variableNum
    if i<=9
        line1=[line1,num2str(i),'. ',variableName{i},'  ;'];
    else
        if i<=18
            line2=[line2,num2str(i),'. ',variableName{i},'  ;'];
        else
            if i<=27
                line3=[line3,num2str(i),'. ',variableName{i},'  ;'];
            end
        end
    end
end
disp(line1)
disp(line2)
if ~isempty(line3)
    disp(line3)
end
end
