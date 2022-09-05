function areaResult=glueMapExplore(glueMap) 
xSize=1024;
ySize=1024;
timeStep=200:200:9600;
% plotAllMap(glueMap);
areaResult=MapCount(glueMap,xSize,ySize,timeStep);
createGluefigure(areaResult.Type0(:,1), [areaResult.Type0(:,2),areaResult.Type1(:,2),areaResult.Type1(:,6),areaResult.Type0(:,6)], [areaResult.Type0(:,4),areaResult.Type1(:,4),areaResult.Type1(:,5),areaResult.Type0(:,5)])
end
function plotAllMap(glueMap)
dirSave=uigetdir();
cd(dirSave);
mkdir(dirSave,'TimeCount');
mkdir(dirSave,'NumberCount');
Type0MapNumerCount=glueMap.Type0MapNumerCount;
Type1MapNumerCount=glueMap.Type1MapNumerCount;
Type0MapTimeCount=glueMap.Type0MapTimeCount;
Type1MapTimeCount=glueMap.Type1MapTimeCount;
dirSave1=strcat(dirSave,'\','TimeCount');
dirSave2=strcat(dirSave,'\','NumberCount');
dirSave3=strcat(dirSave,'\','NumberLabel');
plotGlueMap(Type0MapTimeCount,Type1MapTimeCount,true,dirSave1);
plotGlueMap(Type0MapNumerCount,Type1MapNumerCount,false,dirSave2);
plotGlueMapNumberLab(Type0MapNumerCount,Type1MapNumerCount,dirSave3);
end
function plotGlueMap(Map0,Map1,logScale,dirSave)
endFrame=size(Map1,3);
if logScale==true
    numColor=200;
    colorHot=colormap(hot(numColor));
    Map1=log10(Map1);
    Map0=log10(Map0);
    maxValue=max(max(Map1(:,:,endFrame)));
    Map1=Map1.*fix(numColor/maxValue);
    Map0=Map0.*fix(numColor/maxValue);
end
if logScale==false
    maxValue=max(max(Map1(:,:,endFrame)))/2;
    colorHot=colormap(hot(maxValue));   
end
mkdir(dirSave,'type0Map');
mkdir(dirSave,'type1Map');
for i=1:endFrame
    h=imshow(Map0(:,:,i),colorHot);
    saveFile1=strcat(dirSave,'\','type0Map\',num2str(i));
    saveas(h,saveFile1,'tif'); 
end
for i=1:endFrame
    h=imshow(Map1(:,:,i),colorHot);
    saveFile2=strcat(dirSave,'\','type1Map\',num2str(i));
    saveas(h,saveFile2,'tif'); 
end
end
function plotGlueMapNumberLab(Map0,Map1,dirSave)
endFrame=size(Map1,3);
maxValue=max(max(Map1(:,:,endFrame)));

mkdir(dirSave,'type0Map');
mkdir(dirSave,'type1Map');
for i=1:maxValue
    
    h=imshow(Map0(:,:,endFrame),[0,i]);
    saveFile1=strcat(dirSave,'\','type0Map\',num2str(i));
    saveas(h,saveFile1,'tif');
end
for i=1:maxValue
  
    h=imshow(Map1(:,:,endFrame),[0,i]);
    saveFile2=strcat(dirSave,'\','type1Map\',num2str(i));
    saveas(h,saveFile2,'tif');
end
end
function areaResult=MapCount(glueMap,xSize,ySize,timeStep)
Type0MapNumerCount=glueMap.Type0MapNumerCount;
Type1MapNumerCount=glueMap.Type1MapNumerCount;
Type0MapTimeCount=glueMap.Type0MapTimeCount;
Type1MapTimeCount=glueMap.Type1MapTimeCount;
area0=areaCount(Type0MapTimeCount,xSize,ySize);
area1=areaCount(Type1MapTimeCount,xSize,ySize);
area2=areaCount(Type0MapNumerCount,xSize,ySize);
area3=areaCount(Type1MapNumerCount,xSize,ySize);
areaResult.Type0=[timeStep',area0];
areaResult.Type1=[timeStep' area1];
areaResult.Type00=[timeStep',area2];
areaResult.Type11=[timeStep',area3];
end
function areaResult=areaCount(Map,xSize,ySize)
allSize=(xSize-4)*(ySize-4);
endFrame=size(Map,3);
area=zeros(1,endFrame);
isP=zeros(1,endFrame);
numObj=zeros(1,endFrame);
sizeObj=zeros(1,endFrame);
for i=1:endFrame
    MapTemp=Map(:,:,i);
    [objectInMap,trueOrFalse]=isPocolation(MapTemp,xSize,ySize);
    numObj(i)=objectInMap(1);
    sizeObj(i)=objectInMap(2);
    area(i)=size(MapTemp(MapTemp>=1),1)/allSize;
    isP(i)=trueOrFalse;
end
areaResult=[area',isP',numObj',sizeObj',area'.*isP'];
end
function [objectInMap,trueOrFalse]=isPocolation(Map,xSize,ySize)
mask=im2bw(Map,0);
CC=bwconncomp(mask);
netGlue=regionprops(CC,'area','PixelList');
areaList=[];
searching=true;
for i=1:size(netGlue,1)
    areaList=[areaList,netGlue(i).Area];
    if searching==true
        pixelList=netGlue(i).PixelList;
        inXedgy=pixelList(pixelList(:,1)<=5,:);
        if isempty(inXedgy)
            trueOrFalse=false;
            continue;
        end
        inXedgy=pixelList(pixelList(:,1)>=xSize-5,:);
        if isempty(inXedgy)
            trueOrFalse=false;
            continue;
        end
        inYedgy=pixelList(pixelList(:,2)<=5,:);
        if isempty(inYedgy)
            trueOrFalse=false;
            continue;
        end
        inYedgy=pixelList(pixelList(:,2)>=ySize-5,:);
        if isempty(inYedgy)
            trueOrFalse=false;
            continue;
        end
    end
    trueOrFalse=true;
    searching=false;
end
objectInMap=[size(netGlue,1),mean(areaList)];
end
function createGluefigure(areaResult1, YMatrix1, YMatrix2)
%CREATEFIGURE(AREARESULT1,YMATRIX1,YMATRIX2)
%  AREARESULT1:  vector of x data
%  YMATRIX1:  matrix of y data
%  YMATRIX2:  matrix of y data

%  Auto-generated by MATLAB on 16-Mar-2011 23:54:52

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YGrid','on',...
    'XTickLabel',{'0','100','200','300','400','500','600','700','800',''},...
    'XTick',[0 1200 2400 3600 4800 6000 7200 8400 9600 10000],...
    'XGrid','on',...
    'Position',[0.139503695881731 0.531994981179423 0.760179514255544 0.426740903387699],...
    'FontSize',14,...
    'Color',[1 0.968627452850342 0.921568632125854]);
box(axes1,'on');
hold(axes1,'all');

% Create ylabel
ylabel('Surface Coverage','FontSize',16);

% Create multiple lines using matrix input to plot
plot1 = plot(areaResult1,YMatrix1,'Parent',axes1,'Marker','o');
set(plot1(1),'Marker','square',...
    'DisplayName','areaResult.Type0(:,2) vs areaResult.Type0(:,1)');
set(plot1(2),'Color',[0 0.498039215803146 0],...
    'DisplayName','areaResult.Type1(:,2) vs areaResult.Type1(:,1)');
set(plot1(3),'MarkerFaceColor',[1 0 0],'Color',[1 0 0],...
    'DisplayName','areaResult.Type1(:,6) vs areaResult.Type1(:,1)');
set(plot1(4),'MarkerFaceColor',[0.749019622802734 0 0.749019622802734],...
    'Color',[0.749019622802734 0 0.749019622802734],...
    'DisplayName','areaResult.Type0(:,6) vs areaResult.Type0(:,1)');

% Create axes
axes2 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
    'YMinorGrid','on',...
    'YGrid','on',...
    'XTickLabel',{'0','100','200','300','400','500','600','700','800',''},...
    'XTick',[0 1200 2400 3600 4800 6000 7200 8400 9600 10000],...
    'XGrid','on',...
    'Position',[0.138331573389652 0.0725658720200752 0.76878035902851 0.4075],...
    'FontSize',14,...
    'Color',[0.945098042488098 0.968627452850342 0.949019610881805]);
box(axes2,'on');
hold(axes2,'all');

% Create xlabel
xlabel('Time [Min]','FontSize',16);

% Create multiple lines using matrix input to semilogy
semilogy1 = semilogy(areaResult1,YMatrix2,'Parent',axes2,'Marker','o');
set(semilogy1(1),'DisplayName','Number object in simple');
set(semilogy1(2),'DisplayName','Number object in others');
set(semilogy1(3),'DisplayName','Mean size of others');
set(semilogy1(4),'DisplayName','Mean size of simple');

% Create legend
legend1 = legend(axes2,'show');
set(legend1,...
    'Position',[0.564395452040829 0.341279799247168 0.276663146779303 0.129652864910079]);

% Create textbox
annotation(figure1,'textbox',...
    [0.434002111932417 0.734257214554578 0.125659978880676 0.055207026348808],...
    'String',{'Simple'},...
    'FontSize',20,...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.462449841605066 0.56871141781681 0.125659978880676 0.055207026348808],...
    'String',{'Others'},...
    'FontSize',20,...
    'LineStyle','none');
end