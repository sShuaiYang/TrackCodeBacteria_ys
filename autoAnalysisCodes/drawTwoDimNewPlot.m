function drawTwoDimNewPlot(allData)
clc
close all
disp('please input the draw model')
disp('1.tracking(draw all plot according to tracking result)')
disp('2.mask(draw basic plot only based on mask)')
disp('notice that fpProduce&growthRate analysis only taken on tracking model')
plotChooseParameter=input('modelName=:');
switch plotChooseParameter
    case 1
        plotChooseParameter='tracking';
        dataAll=allData.trackingResult;
        time=dataAll(:,1);
        maxLene=dataAll(:,2);
        minLen=dataAll(:,3);
        CyOFP=dataAll(:,4);
        GFP=dataAll(:,5);
        mScalet=dataAll(:,6);
        RFP=dataAll(:,7);
        angle=dataAll(:,8);
        area=dataAll(:,9);
        growthRate=dataAll(:,10);
        fpProduceGFP=dataAll(:,11);
        treeSize=dataAll(:,12);
        generation=dataAll(:,13);
        CyPet=dataAll(:,14);
        Venus=dataAll(:,15);
        mAmetrine=dataAll(:,16);
        fpProduceCyOFP=dataAll(:,17);
        fpProducemScalet=dataAll(:,18);
        fpProduceRFP=dataAll(:,19);
        fpProduceCyPete=dataAll(:,20);
        fpProduceVenus=dataAll(:,21);
        fpProducemAmetrine=dataAll(:,22);
        redControl=dataAll(:,23);
        blueControl=dataAll(:,24);
        greenControl=dataAll(:,25);
    case 2
        plotChooseParameter='mask';
        dataAll=allData.maskResult;
        time=dataAll(:,1);
        maxLene=dataAll(:,2);
        minLen=dataAll(:,3);
        CyOFP=dataAll(:,4);
        GFP=dataAll(:,5);
        mScalet=dataAll(:,6);
        RFP=dataAll(:,7);
        CyPet=dataAll(:,8);
        Venus=dataAll(:,9);
        mAmetrine=dataAll(:,10);
        redControl=dataAll(:,11);
        blueControl=dataAll(:,12);
        greenControl=dataAll(:,13);
end
switch plotChooseParameter
    case 1
        plotChooseParameter='tracking';
    case 2
        plotChooseParameter='mask';
end
variableName=checkData(allData,plotChooseParameter);
variValue1=input('please type the variable____x:');
iLine=[];
jLine=[];
if isnumeric(variValue1);
    variValue1=variableName{variValue1};
else
    iLine=eval(variValue1);
end
variValue2=input('please type the variable____y:');
if isnumeric(variValue2);
    variValue2=variableName{variValue2};
else
    jLine=eval(variValue2);
end
plotStyle=input('please input the plotStyle    1.dot   2.gaussian:');
needPara=input('need seperate according to other parameter?  [y/n] :');
if strcmp(needPara,'y')
    otherParameter=input('please type the other parameter:');
    otherParameter=variableName{otherParameter};
    paraRange=input('please input the divide detail [  diff/  n  piece /specify region ]:');
end
switch plotChooseParameter
    case 'tracking'
        if isempty(iLine)
            iLine=getTrackingResultLine(allData.trackingResult,variValue1);
        end
        if isempty(jLine)
            jLine=getTrackingResultLine(allData.trackingResult,variValue2);
        end
        iLine(isnan(jLine))=[];
        jLine(isnan(jLine))=[];
        jLine(isnan(iLine))=[];
        iLine(isnan(iLine))=[];
        if strcmp(needPara,'y');
            if strcmp(otherParameter,'tag')
                cLine=allData.trackingResultTag;
                tagA=unique(cLine);
                xNum=ceil(sqrt(numel(tagA)));
                yNum=ceil(numel(tagA)/xNum);
                for i=1:numel(tagA)
                    iLineNew=iLine(strcmp(cLine,tagA{i}));
                    jLineNew=jLine(strcmp(cLine,tagA{i}));
                    switch plotStyle
                        case 1
                            hold on
                            subplot(xNum,yNum,i)
                            plot(iLineNew,jLineNew,'DisplayName',tagA{i},'linestyle','none','marker','o','markerSize',6)
                            set(gca,'FontSize',20);
                            xlabel(variValue1,'FontSize',20);
                            ylabel(variValue2,'FontSize',20);
                            legend(gca,'show');
                        case 2
                            hold on
                            subplot(xNum,yNum,i)
                            result=twoDimGaussian(iLineNew,jLineNew);
                            set(gca,'FontSize',20);
                            xlabel(variValue1,'FontSize',20);
                            ylabel(variValue2,'FontSize',20);
                            legend(gca,'show');
                    end
                end
            else
                switch otherParameter
                    case 'tagValue'
                        cLine=allData.trackingResultTagValue;
                    otherwise
                        cLine=getTrackingResultLine(allData.trackingResult,otherParameter);
                end
                if strcmp(paraRange,'diff')
                    tagA=unique(cLine);
                    tagA(isnan(tagA))=[];
                    xNum=ceil(sqrt(numel(tagA)));
                    yNum=ceil(numel(tagA)/xNum);
                    for i=1:numel(tagA)
                        iLineNew=iLine(cLine==tagA{i});
                        jLineNew=jLine(cLine==tagA{i});
                        switch plotStyle
                            case 1
                                hold on
                                subplot(xNum,yNum,i)
                                plot(iLineNew,jLineNew,'DisplayName',num2str(tagA(i)),'linestyle','none','marker','o','markerSize',6)
                                set(gca,'FontSize',20);
                                xlabel(variValue1,'FontSize',20);
                                ylabel(variValue2,'FontSize',20);
                                legend(gca,'show');
                            case 2
                                hold on
                                subplot(xNum,yNum,i)
                                result=twoDimGaussian(iLineNew,jLineNew);
                                set(gca,'FontSize',20);
                                xlabel(variValue1,'FontSize',20);
                                ylabel(variValue2,'FontSize',20);
                                legend(gca,'show');
                        end
                    end
                else
                    if isnumeric(paraRange) && numel(paraRange)==1
                        maxValue=max(cLine);
                        minValue=min(cLine);
                        divRange=linspace(minValue,maxValue,paraRange+1);
                        xNum=ceil(sqrt(numel(divRange)-1));
                        yNum=ceil((numel(divRange)-1)/xNum);
                        for i=1:numel(divRange)-1
                            iLineNew=iLine(cLine>=divRange(i) & cLine<divRange(i+1));
                            jLineNew=jLine(cLine>=divRange(i) & cLine<divRange(i+1));
                            switch plotStyle
                                case 1
                                    hold on
                                    subplot(xNum,yNum,i)
                                    plot(iLineNew,jLineNew,'DisplayName',[otherParameter,':[',num2str(divRange(i)),'--',num2str(divRange(i+1)),']'],'linestyle','none','marker','o','markerSize',6)
                                    set(gca,'FontSize',20);
                                    xlabel(variValue1,'FontSize',20);
                                    ylabel(variValue2,'FontSize',20);
                                    legend(gca,'show');
                                case 2
                                    hold on
                                    subplot(xNum,yNum,i)
                                    result=twoDimGaussian(iLineNew,jLineNew);
                                    set(gca,'FontSize',20);
                                    xlabel(variValue1,'FontSize',20);
                                    ylabel(variValue2,'FontSize',20);
                                    legend(gca,'show');
                            end
                        end
                    else
                        divRange=paraRange;
                        xNum=ceil(sqrt(numel(divRange)-1));
                        yNum=ceil((numel(divRange)-1)/xNum);
                        for i=1:numel(divRange)-1
                            iLineNew=iLine(cLine>=divRange(i) & cLine<divRange(i+1));
                            jLineNew=jLine(cLine>=divRange(i) & cLine<divRange(i+1));
                            switch plotStyle
                                case 1
                                    hold on
                                    subplot(xNum,yNum,i)
                                    plot(iLineNew,jLineNew,'DisplayName',[otherParameter,':[',num2str(divRange(i)),'--',num2str(divRange(i+1)),']'],'linestyle','none','marker','o','markerSize',6)
                                    set(gca,'FontSize',20);
                                    xlabel(variValue1,'FontSize',20);
                                    ylabel(variValue2,'FontSize',20);
                                    legend(gca,'show');
                                case 2
                                    hold on
                                    subplot(xNum,yNum,i)
                                    result=twoDimGaussian(iLineNew,jLineNew);
                                    set(gca,'FontSize',20);
                                    xlabel(variValue1,'FontSize',20);
                                    ylabel(variValue2,'FontSize',20);
                                    legend(gca,'show');
                            end
                        end
                    end
                end
            end
        else
            switch plotStyle
                case 1
                    hold on
                    plot(iLine,jLine,'linestyle','none','marker','o','markerSize',6)
                    set(gca,'FontSize',20);
                    xlabel(variValue1,'FontSize',20);
                    ylabel(variValue2,'FontSize',20);
                case 2
                    hold on
                    result=twoDimGaussian(iLine,jLine);
                    set(gca,'FontSize',20);
                    xlabel(variValue1,'FontSize',20);
                    ylabel(variValue2,'FontSize',20);
            end
        end
    case 'mask'
        if isempty(iLine)
            iLine=getMaskResultLine(allData.maskResult,variValue1);
        end
        if isempty(jLine)
            jLine=getMaskResultLine(allData.maskResult,variValue2);
        end
        iLine(isnan(jLine))=[];
        jLine(isnan(jLine))=[];
        jLine(isnan(iLine))=[];
        iLine(isnan(iLine))=[];
        if strcmp(needPara,'y');
            if strcmp(otherParameter,'tag')
                cLine=allData.maskResultTag;
                tagA=unique(cLine);
                xNum=ceil(sqrt(numel(tagA)));
                yNum=ceil(numel(tagA)/xNum);
                for i=1:numel(tagA)
                    iLineNew=iLine(strcmp(cLine,tagA{i}));
                    jLineNew=jLine(strcmp(cLine,tagA{i}));
                    switch plotStyle
                        case 1
                            hold on
                            subplot(xNum,yNum,i)
                            plot(iLineNew,jLineNew,'DisplayName',tagA{i},'linestyle','none','marker','o','markerSize',6)
                            set(gca,'FontSize',20);
                            xlabel(variValue1,'FontSize',20);
                            ylabel(variValue2,'FontSize',20);
                            legend(gca,'show');
                        case 2
                            hold on
                            subplot(xNum,yNum,i)
                            result=twoDimGaussian(iLineNew,jLineNew);
                            set(gca,'FontSize',20);
                            xlabel(variValue1,'FontSize',20);
                            ylabel(variValue2,'FontSize',20);
                            legend(gca,'show');
                    end
                end
            else
                switch otherParameter
                    case 'tagValue'
                        cLine=allData.maskResultTagValue;
                    otherwise
                        cLine=getMaskResultLine(allData.maskResult,otherParameter);
                end
                if strcmp(paraRange,'diff')
                    tagA=unique(cLine);
                    tagA(isnan(tagA))=[];
                    xNum=ceil(sqrt(numel(tagA)));
                    yNum=ceil(numel(tagA)/xNum);
                    for i=1:numel(tagA)
                        iLineNew=iLine(cLine==tagA{i});
                        jLineNew=jLine(cLine==tagA{i});
                        switch plotStyle
                            case 1
                                hold on
                                subplot(xNum,yNum,i)
                                plot(iLineNew,jLineNew,'DisplayName',num2str(tagA(i)),'linestyle','none','marker','o','markerSize',6)
                                set(gca,'FontSize',20);
                                xlabel(variValue1,'FontSize',20);
                                ylabel(variValue2,'FontSize',20);
                                legend(gca,'show');
                            case 2
                                hold on
                                subplot(xNum,yNum,i)
                                result=twoDimGaussian(iLineNew,jLineNew);
                                set(gca,'FontSize',20);
                                xlabel(variValue1,'FontSize',20);
                                ylabel(variValue2,'FontSize',20);
                                legend(gca,'show');
                        end
                    end
                else
                    if isnumeric(paraRange) && numel(paraRange)==1
                        maxValue=max(cLine);
                        minValue=min(cLine);
                        divRange=linspace(minValue,maxValue,paraRange+1);
                        xNum=ceil(sqrt(numel(divRange)-1));
                        yNum=ceil((numel(divRange)-1)/xNum);
                        for i=1:numel(divRange)-1
                            iLineNew=iLine(cLine>=divRange(i) & cLine<divRange(i+1));
                            jLineNew=jLine(cLine>=divRange(i) & cLine<divRange(i+1));
                            switch plotStyle
                                case 1
                                    hold on
                                    subplot(xNum,yNum,i)
                                    plot(iLineNew,jLineNew,'DisplayName',[otherParameter,':[',num2str(divRange(i)),'--',num2str(divRange(i+1)),']'],'linestyle','none','marker','o','markerSize',6)
                                    set(gca,'FontSize',20);
                                    xlabel(variValue1,'FontSize',20);
                                    ylabel(variValue2,'FontSize',20);
                                    legend(gca,'show');
                                case 2
                                    hold on
                                    subplot(xNum,yNum,i)
                                    result=twoDimGaussian(iLineNew,jLineNew);
                                    set(gca,'FontSize',20);
                                    xlabel(variValue1,'FontSize',20);
                                    ylabel(variValue2,'FontSize',20);
                                    legend(gca,'show');
                            end
                        end
                    else
                        divRange=paraRange;
                        xNum=ceil(sqrt(numel(divRange)-1));
                        yNum=ceil((numel(divRange)-1)/xNum);
                        for i=1:numel(divRange)-1
                            iLineNew=iLine(cLine>=divRange(i) & cLine<divRange(i+1));
                            jLineNew=jLine(cLine>=divRange(i) & cLine<divRange(i+1));
                            switch plotStyle
                                case 1
                                    hold on
                                    subplot(xNum,yNum,i)
                                    plot(iLineNew,jLineNew,'DisplayName',[otherParameter,':[',num2str(divRange(i)),'--',num2str(divRange(i+1)),']'],'linestyle','none','marker','o','markerSize',6)
                                    set(gca,'FontSize',20);
                                    xlabel(variValue1,'FontSize',20);
                                    ylabel(variValue2,'FontSize',20);
                                    legend(gca,'show');
                                case 2
                                    hold on
                                    subplot(xNum,yNum,i)
                                    result=twoDimGaussian(iLineNew,jLineNew);
                                    set(gca,'FontSize',20);
                                    xlabel(variValue1,'FontSize',20);
                                    ylabel(variValue2,'FontSize',20);
                                    legend(gca,'show');
                            end
                        end
                    end
                end
            end
        else
            switch plotStyle
                case 1
                    hold on
                    plot(iLine,jLine,'linestyle','none','marker','o','markerSize',6)
                    set(gca,'FontSize',20);
                    xlabel(variValue1,'FontSize',20);
                    ylabel(variValue2,'FontSize',20);
                case 2
                    hold on
                    result=twoDimGaussian(iLine,jLine);
                    set(gca,'FontSize',20);
                    xlabel(variValue1,'FontSize',20);
                    ylabel(variValue2,'FontSize',20);
            end
        end
end
end
function iLine=getMaskResultLine(dataAll,inputVari)
switch inputVari
    case 'time'
        iLine=dataAll(:,1);
    case 'maxLen'
        iLine=dataAll(:,2);
    case 'minLen'
        iLine=dataAll(:,3);
    case 'CyOFP'
        iLine=dataAll(:,4);
    case 'GFP'
        iLine=dataAll(:,5);
    case 'mScalet'
        iLine=dataAll(:,6);
    case 'RFP'
        iLine=dataAll(:,7);
    case 'CyPet'
        iLine=dataAll(:,8);
    case 'Venus'
        iLine=dataAll(:,9);
    case 'mAmetrine'
        iLine=dataAll(:,10);
    case 'redControl'
        iLine=dataAll(:,11);
    case 'blueControl'
        iLine=dataAll(:,12);
    case 'greenControl'
        iLine=dataAll(:,13);
    otherwise
        iLine=[];
end
end
function iLine=getTrackingResultLine(dataAll,inputVari)
switch inputVari
    case 'time'
        iLine=dataAll(:,1);
    case 'maxLen'
        iLine=dataAll(:,2);
    case 'minLen'
        iLine=dataAll(:,3);
    case 'CyOFP'
        iLine=dataAll(:,4);
    case 'GFP'
        iLine=dataAll(:,5);
    case 'mScalet'
        iLine=dataAll(:,6);
    case 'RFP'
        iLine=dataAll(:,7);
    case 'angle'
        iLine=dataAll(:,8);
    case 'area'
        iLine=dataAll(:,9);
    case 'growthRate'
        iLine=dataAll(:,10);
    case 'fpProduceGFP'
        iLine=dataAll(:,11);
    case 'treeSize'
        iLine=dataAll(:,12);
    case 'generation'
        iLine=dataAll(:,13);
    case 'CyPet'
        iLine=dataAll(:,14);
    case 'Venus'
        iLine=dataAll(:,15);
    case 'mAmetrine'
        iLine=dataAll(:,16);
    case 'fpProduceCyOFP'
        iLine=dataAll(:,17);
    case 'fpProducemScalet'
        iLine=dataAll(:,18);
    case 'fpProduceRFP'
        iLine=dataAll(:,19);
    case 'fpProduceCyPet'
        iLine=dataAll(:,20);
    case 'fpProduceVenus'
        iLine=dataAll(:,21);
    case 'fpProducemAmetrine'
        iLine=dataAll(:,22);
    case 'redControl'
        iLine=dataAll(:,23);
    case 'blueControl'
        iLine=dataAll(:,24);
    case 'greenControl'
        iLine=dataAll(:,25);
end
end
function variableName=checkData(allData,plotChooseParameter)
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
variableName{16}='treeSize';
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
maskImage=allData.maskResult;
trackingResult=allData.trackingResult;
switch plotChooseParameter
    case 'tracking'
        chooseNum=27;
    case 'mask'
        variableName=variableName([1:10,14:15,25:27]);
        chooseNum=15;
        chooseIndex=true(chooseNum,1);
end
chooseIndex=true(chooseNum,1);
for i=4:10
    if all(isnan(maskImage(:,i))) 
        chooseIndex(i)=0;
        chooseIndex(i+14)=0;
    end
end
if ~isempty(trackingResult)
    for i=25:27
        if all(isnan(trackingResult(:,i-2)))
            chooseIndex(i)=0;
        end
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
disp('now please input the variables');
disp(line1)
disp(line2)
if ~isempty(line3)
    disp(line3)
end
end
function result=twoDimGaussian(a,b)
% 二维高斯blur
% xMin=min(a);
% yMin=min(b);
% xMax=max(a);
% yMax=max(b);
xMin=0;
xMax=60;
yMin=0;
yMax=400;
sigma1=(xMax-xMin)/50;
interval1=(xMax-xMin)/200;
sigma2=(yMax-yMin)/50;
interval2=(yMax-yMin)/200;
rowPos=xMin:interval1:xMax;
rowPos=rowPos';
columnPos=yMin:interval2:yMax;
resultImage=zeros(numel(rowPos),numel(columnPos));
for i=1:numel(a)
%     disp(i)
    rowGauss=1/(2*pi)^0.5/sigma1*exp(-(rowPos-a(i)).^2/2/sigma1^2);
    columnGauss=1/(2*pi)^0.5/sigma2*exp(-(columnPos-b(i)).^2/2/sigma2^2);
    resultImage=resultImage+rowGauss*columnGauss;
end
result.Image=resultImage;
result.numPoint=numel(a);
result.x=rowPos';
result.y=columnPos;
surface(result.x,result.y,result.Image','LineStyle','none');
xlim(gca,[xMin,xMax])
ylim(gca,[yMin,yMax])
end