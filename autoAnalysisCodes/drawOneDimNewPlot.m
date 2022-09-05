function drawOneDimNewPlot(allData)
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
variableName=checkData(allData,plotChooseParameter);
variValue=input('please type the variable____:');
iLine=[];
if isnumeric(variValue);
    variValue=variableName{variValue};
else
    iLine=eval(variValue);
end
mainVara=variValue;
needPara=input('need seperate according to other parameter?  [y/n] :');
if strcmp(needPara,'y')
    otherParameter=input('please type the other parameter:');
    otherParameter=variableName{otherParameter};
    paraRange=input('please input the divide detail [  diff/  n  piece /specify region ]:');
end
switch plotChooseParameter
    case 'tracking'
        if isempty(iLine)
            iLine=getTrackingResultLine(allData.trackingResult,mainVara);
        end
        if strcmp(needPara,'y');
            if strcmp(otherParameter,'tag')
                cLine=allData.trackingResultTag;
                tagA=unique(cLine);
                colorAll=colormap(jet(numel(tagA)));
                divNum=input('please input hist div number::');
                for i=1:numel(tagA)
                    kLine=iLine(strcmp(cLine,tagA{i}));
                    kLine=kLine(~isnan(kLine));
                    if ~isempty(kLine)
                        [b,a]=hist(kLine,divNum);
                        result{i}=[a;b];
                        hold on;plot(a,b./sum(b),'DisplayName',tagA{i},'color',colorAll(i,:),'linestyle','none','marker','o','markerSize',6)
                    end
                end
                legend(gca,'show');
            else
                switch otherParameter
                    case 'tagValue'
                        cLine=allData.trackingResultTagValue;
                    otherwise
                        cLine=getTrackingResultLine(allData.trackingResult,otherParameter);
                end
                if strcmp(paraRange,'diff')
                    tagA=unique(cLine);tagA(isnan(tagA))=[];
                    colorAll=colormap(jet(numel(tagA)));
                    divNum=input('please input hist div number:');
                    for i=1:numel(tagA)
                        kLine=iLine(cLine==tagA(i));
                        kLine=kLine(~isnan(kLine));
                        if ~isempty(kLine)
                            [b,a]=hist(kLine,divNum);
                            result{i}=[a;b];
                            hold on;plot(a,b./sum(b),'DisplayName',num2str(tagA(i)),'color',colorAll(i,:),'linestyle','none','marker','o','markerSize',6)
                        end
                    end
                    legend(gca,'show');
                else
                    if isnumeric(paraRange) && numel(paraRange)==1
                        maxValue=max(cLine);
                        minValue=min(cLine);
                        divRange=linspace(minValue,maxValue,paraRange+1);
                        colorAll=colormap(jet(numel(divRange)));
                        divNum=input('please input hist div number:');
                        for i=1:numel(divRange)-1
                            kLine=iLine(cLine>=divRange(i) & cLine<divRange(i+1));
                            kLine=kLine(~isnan(kLine));
                            if ~isempty(kLine)
                                [b,a]=hist(kLine,divNum);
                                result{i}=[a;b];
                                hold on;plot(a,b./sum(b),'DisplayName',[otherParameter,':[',num2str(divRange(i)),'--',num2str(divRange(i+1)),']'],'color',colorAll(i,:),'linestyle','none','marker','o','markerSize',6);
                            end
                        end
                        legend(gca,'show');
                    else
                        divRange=paraRange;
                        colorAll=colormap(jet(numel(divRange)-1));
                        divNum=input('please input hist div number:');
                        for i=1:numel(divRange)-1
                            kLine=iLine(cLine>=divRange(i) & cLine<divRange(i+1));
                            kLine=kLine(~isnan(kLine));
                            if ~isempty(kLine)
                                [b,a]=hist(kLine,divNum);
                                result{i}=[a;b];
                                hold on;plot(a,b./sum(b),'DisplayName',[otherParameter,':[',num2str(divRange(i)),'--',num2str(divRange(i+1)),']'],'color',colorAll(i,:),'linestyle','none','marker','o','markerSize',6);
                            end
                        end
                        legend(gca,'show');
                    end
                end
            end
        else
            divNum=input('please input hist div number:');
            [b,a]=hist(iLine,divNum);
            plot(a,b./sum(b),'linestyle','none','marker','o','markerSize',6)
            result{1}=[a;b];
        end
    case 'mask'
        if isempty(iLine)
            iLine=getMaskResultLine(allData.maskResult,mainVara);
        end
        if strcmp(needPara,'y');
            if strcmp(otherParameter,'tag')
                cLine=allData.maskResultTag;
                cLine=cLine(iLine<=100 & iLine>=0);
                iLine=iLine(iLine<=100 & iLine>=0);
                tagA=unique(cLine);
                colorAll=colormap(jet(numel(tagA)));
                divNum=input('please input hist div number:');
                for i=1:numel(tagA)
                    kLine=iLine(strcmp(cLine,tagA{i}));
                    kLine=kLine(~isnan(kLine));
                    if ~isempty(kLine)
                        [b,a]=hist(kLine,divNum);
                        result{i}=[a;b];
                        hold on;plot(a,b./sum(b),'DisplayName',tagA{i},'color',colorAll(i,:),'linestyle','none','marker','o','markerSize',6)
                    end
                end
                legend(gca,'show');
            else
                switch otherParameter
                    case 'tagValue'
                        cLine=allData.maskResultTagValue;
                    otherwise
                        cLine=getMaskResultLine(allData.maskResult,otherParameter);
                end
                if strcmp(paraRange,'diff')
                    tagA=unique(cLine);tagA(isnan(tagA))=[];
                    colorAll=colormap(jet(numel(tagA)));
                    divNum=input('please input hist div number:');
                    for i=1:numel(tagA)
                        kLine=iLine(cLine==tagA(i));
                        kLine=kLine(~isnan(kLine));
                        if ~isempty(kLine)
                            [b,a]=hist(kLine,divNum);
                            result{i}=[a;b];
                            hold on;plot(a,b./sum(b),'DisplayName',num2str(tagA(i)),'color',colorAll(i,:),'linestyle','none','marker','o','markerSize',6)
                        end
                    end
                    legend(gca,'show');
                end
                if isnumeric(paraRange) && numel(paraRange)==1
                    maxValue=max(cLine);
                    minValue=min(cLine);
                    divRange=linspace(minValue,maxValue,paraRange+1);
                    colorAll=colormap(jet(numel(divRange)));
                    divNum=input('please input hist div number:');
                    for i=1:numel(divRange)-1
                        kLine=iLine(cLine>=divRange(i) & cLine<divRange(i+1));
                        kLine=kLine(~isnan(kLine));
                        if ~isempty(kLine)
                            [b,a]=hist(kLine,divNum);
                            result{i}=[a;b];
                            hold on;plot(a,b./sum(b),'DisplayName',[otherParameter,':[',num2str(divRange(i)),'--',num2str(divRange(i+1)),']'],'color',colorAll(i,:),'linestyle','none','marker','o','markerSize',6);
                        end
                    end
                    legend(gca,'show');
                else
                    divRange=paraRange;
                    colorAll=colormap(jet(numel(divRange)-1));
                    divNum=input('please input hist div number:');
                    for i=1:numel(divRange)-1
                        kLine=iLine(cLine>=divRange(i) & cLine<divRange(i+1));
                        kLine=kLine(~isnan(kLine));
                        if ~isempty(kLine)
                            [b,a]=hist(kLine,divNum);
                            result{i}=[a;b];
                            hold on;plot(a,b./sum(b),'DisplayName',[otherParameter,':[',num2str(divRange(i)),'--',num2str(divRange(i+1)),']'],'color',colorAll(i,:),'linestyle','none','marker','o','markerSize',6);
                        end
                    end
                    legend(gca,'show');
                end
            end
        else
            divNum=input('please input hist div number:');
            [b,a]=hist(iLine,divNum);
            plot(a,b./sum(b),'linestyle','none','marker','o','markerSize',6)
            result{1}=[a;b];
        end
end
xlabel(variValue,'FontSize',20);
needFit=input('do you want to fit your data?  [y/n] :');
switch needFit
    case 'y'
        fitType=input('input your fit type gauss/gamma?  :');
        switch fitType
            case 'gauss'
                for i=1:numel(result)
                    try
                        aveNum(i)=prepareGaussCurveData(result{i}(1,:),result{i}(2,:));
                    end
                end
            case 'gamma'
                for i=1:numel(result)
                    try
                        prepareGammaCurveData(result{i}(1,:),result{i}(2,:));
                    end
                end
        end
end
end
function aveNum=prepareGaussCurveData(x,y)
y=y./sum(y);
% Set up fittype and options.
ft = fittype('gauss1');
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf];
opts.StartPoint = [max(y),mean(x),mean(x)/10];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( x',y', ft, opts );

% Plot fit with data.
hold on;
h = plot( fitresult);
text(fitresult.b1,fitresult.a1,['a=',num2str(fitresult.b1),'|b=',num2str(fitresult.c1)])
aveNum=fitresult.b1;
end
function prepareGammaCurveData(x,y)
y=y./sum(y);
meanX=x(y==max(y));
x1=meanX(1)/4;
x=x/x1;
y(x>=20)=[];
x(x>=20)=[];
% Set up fittype and options.
ft = fittype( 'b^(-a)*c/gamma(a)*x^(a-1)*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [0 0 0];
opts.StartPoint = [10,0.4,max(y)];
opts.Upper = [100 100 10];

% Fit model to data.
[fitresult, gof] = fit(x',y',ft,opts);
hold on;
x=linspace(min(x),max(x),1000);
h = plot(x*x1,fitresult.b^(-fitresult.a)*fitresult.c/gamma(fitresult.a)*x.^(fitresult.a-1).*exp(-x/fitresult.b));
text(fitresult.a*fitresult.b*x1,fitresult.b^(-fitresult.a)*fitresult.c/gamma(fitresult.a)*(fitresult.a*fitresult.b).^(fitresult.a-1).*exp(-(fitresult.a*fitresult.b/fitresult.b)),['k1/γ2=',num2str(fitresult.a),'|k2/γ1=',num2str(fitresult.b*x1/log(2))]);
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
        chooseIndex=true(chooseNum,1);
        for i=25:27
            if all(isnan(trackingResult(:,i-2)))
                chooseIndex(i)=0;
            end
        end
    case 'mask'
        variableName=variableName([1:10,14:15,25:27]);
        chooseNum=15;
        chooseIndex=true(chooseNum,1);
end
for i=4:10
    if all(isnan(maskImage(:,i))) 
        chooseIndex(i)=0;
        chooseIndex(i+14)=0;
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
