function [ bioTree ] = bioTreeAllInfoGet( bioTree,dirFile)
%BIOTREEALLINFOGET Summary of this function goes here
%   Detailed explanation goes here
load([dirFile,'\Tracking\frameInfo.mat'])
bioTreeInitialTime=frameInfo(1,1:6);
frameInfoNew(1:2:numel(bioTree)-1,:)=frameInfo(1:numel(bioTree)/2,:);
frameInfoNew(2:2:numel(bioTree),:)=frameInfo(1:numel(bioTree)/2,:);
% frameInfoNew=frameInfo;
% bestPositionAccumulation=bioTree{1}.bestPositionAccumulation;
for i=1:numel(bioTree)
    bioTreeTimer(i)=etime(frameInfoNew(i,1:6),bioTreeInitialTime)/60;   %min
end
bioTree{1}.bioTreeInitialTime=bioTreeInitialTime;
bioTree{1}.bioTreeTimer=bioTreeTimer;
fluoChannel{1,1}='CyOFP';
fluoChannel{2,1}='sfGFP';
fluoChannel{3,1}='mScarletI';
fluoChannel{4,1}='TDsmURFP';
fluoChannel{5,1}='CyPet';
fluoChannel{6,1}='Venus';
fluoChannel{7,1}='mAmetrine';
bioTree=bioTreeAllChannel(bioTree);
for iChannel=1:numel(fluoChannel)
    dirFluoFile=[dirFile,'\',fluoChannel{iChannel},'new'];
    try
        load([dirFluoFile,'\frameInfo.mat'])
    catch err
        continue
    end
    bioTreeFrame=[];
    nameList=dir(dirFluoFile);
    for i=1:size(frameInfo,1)
        bioTreeTimer=etime(frameInfo(i,1:6),bioTree{1}.bioTreeInitialTime)/60;
        diffTime=abs(bioTree{1}.bioTreeTimer-bioTreeTimer);
        [~,index]=find(diffTime==min(diffTime));
        bioTreeFrame(i)=index(2);
%         bioTreeFrame(i)=index(1);
    end
    [~,checkFrame]=find(diff(bioTreeFrame)==0);
    bioTreeFrame(checkFrame+1)=0;
    for iPic=1:numel(bioTreeFrame)
        if bioTreeFrame(iPic)==0
            continue
        end
        %         image=load([dirFluoFile,'\image',fluoChannel{iChannel},num2str(iPic,'%05.f'),'.mat']);
        image=load([dirFluoFile,'\image','Tacking',num2str(iPic,'%05.f'),'.mat']);
        image=image.imageTracking;
%         image=image.imageI;
%         se=translate(strel(1),bestPositionAccumulation(bioTreeFrame(iPic),:));
%         image=imdilate(image,se);
        for iframe=1:bioTreeFrame(iPic)
            for iRoot=1:numel(bioTree{iframe}.root)
                traceInfo=bioTree{iframe}.root{iRoot}.traceInfo;
                bioTree{iframe}.root{iRoot}.linkTimer=bioTree{1}.bioTreeTimer(iframe:iframe+numel(traceInfo.pixelIdxList)-1)';
                if numel(traceInfo.pixelIdxList)+iframe>bioTreeFrame(iPic)
                    tracePixel=traceInfo.pixelIdxList{bioTreeFrame(iPic)-iframe+1};
                    switch iChannel
                        case 1
                            bioTree{iframe}.root{iRoot}.CyOFP(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                        case 2
                            bioTree{iframe}.root{iRoot}.GFP(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                        case 3
                            bioTree{iframe}.root{iRoot}.mScalet(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                        case 4
                            bioTree{iframe}.root{iRoot}.RFP(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                        case 5
                            bioTree{iframe}.root{iRoot}.CyPet(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                        case 6
                            bioTree{iframe}.root{iRoot}.Venus(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                        case 7
                            bioTree{iframe}.root{iRoot}.mAmetrine(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                    end
                end
            end
            for iNode=1:numel(bioTree{iframe}.node)
                for iOut=1:numel(bioTree{iframe}.node{iNode}.Out)
                    traceInfo=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo;
                    bioTree{iframe}.node{iNode}.Out{iOut}.linkTimer=bioTree{1}.bioTreeTimer(iframe:iframe+numel(traceInfo.pixelIdxList)-1)';
                    if numel(traceInfo.pixelIdxList)+iframe>bioTreeFrame(iPic)
                        tracePixel=traceInfo.pixelIdxList{bioTreeFrame(iPic)-iframe+1};
                        switch iChannel
                            case 1
                                bioTree{iframe}.node{iNode}.Out{iOut}.CyOFP(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                            case 2
                                bioTree{iframe}.node{iNode}.Out{iOut}.GFP(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                            case 3
                                bioTree{iframe}.node{iNode}.Out{iOut}.mScalet(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                            case 4
                                bioTree{iframe}.node{iNode}.Out{iOut}.RFP(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                            case 5
                                bioTree{iframe}.node{iNode}.Out{iOut}.CyPet(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                            case 6
                                bioTree{iframe}.node{iNode}.Out{iOut}.Venus(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                            case 7
                                bioTree{iframe}.node{iNode}.Out{iOut}.mAmetrine(bioTreeFrame(iPic)-iframe+1)=mean(image(tracePixel));
                        end
                    end
                end
            end
        end
    end
end

fluoChannel{8,1}='RedMask';
fluoChannel{9,1}='BlueMask';
fluoChannel{10,1}='GreenMask';
% mip=bioTree{1}.mip;
% for iChannel=8:10
%     if isempty(mip.Calib.lightControl);
%         return
%     end
%     if strcmp(mip.Calib.lightControl,'static')
%         dirFluoFile=[dirFile,'\',fluoChannel{iChannel}];
%         try
%             load([dirFluoFile,'\frameInfo.mat'])
%         catch err
%             return
%         end
%         staticPower=frameInfo(1,15);
%         staticLightBackGround=frameInfo(1,14);
%         staticBackGroundImage=import_tiff_stack([dirFile,'\',fluoChannel{iChannel}(1:end-4),'Control\image',fluoChannel{iChannel}(1:end-4),'Control00001.tif']);
%         for iframe=1:numel(bioTree)
%             for iRoot=1:numel(bioTree{iframe}.root)
%                 traceInfo=bioTree{iframe}.root{iRoot}.traceInfo;
%                 for iTrace=1:numel(traceInfo.pixelIdxList);
%                     tracePixel=traceInfo.pixelIdxList{iTrace};
%                     switch iChannel
%                         case 8
%                             bioTree{iframe}.root{iRoot}.RedControl(iTrace)=mean(double(staticBackGroundImage(tracePixel)))*16/staticLightBackGround*staticPower;
%                         case 9
%                             bioTree{iframe}.root{iRoot}.BlueControl(iTrace)=mean(double(staticBackGroundImage(tracePixel)))*16/staticLightBackGround*staticPower;
%                         case 10
%                             bioTree{iframe}.root{iRoot}.GreenControl(iTrace)=mean(double(staticBackGroundImage(tracePixel)))*16/staticLightBackGround*staticPower;
%                             
%                     end
%                 end
%             end
%             for iNode=1:numel(bioTree{iframe}.node)
%                 for iOut=1:numel(bioTree{iframe}.node{iNode}.Out)
%                     traceInfo=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo;
%                     for iTrace=1:numel(traceInfo.pixelIdxList)
%                         tracePixel=traceInfo.pixelIdxList{iTrace};
%                         switch iChannel
%                             case 8
%                                 bioTree{iframe}.node{iNode}.Out{iOut}.RedControl(iTrace)=mean(double(staticBackGroundImage(tracePixel)))*16/staticLightBackGround*staticPower;
%                             case 9
%                                 bioTree{iframe}.node{iNode}.Out{iOut}.BlueControl(iTrace)=mean(double(staticBackGroundImage(tracePixel)))*16/staticLightBackGround*staticPower;
%                             case 10
%                                 bioTree{iframe}.node{iNode}.Out{iOut}.GreenControl(iTrace)=mean(double(staticBackGroundImage(tracePixel)))*16/staticLightBackGround*staticPower;
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     if strcmp(mip.Calib.lightControl,'static')
%         dirFluoFile=[dirFile,'\',fluoChannel{iChannel}];
%         try
%             load([dirFluoFile,'\frameInfo.mat'])
%         catch err
%             return
%         end
%         for i=1:size(frameInfo,1)-1
%             deltaTime(i)=abs(etime(frameInfo(i+1,1:6),frameInfo(i,1:6)));
%         end
%         deltaTime(end+1)=deltaTime(end);
%         frameInfo(:,15)=frameInfo(:,15).*frameInfo(:,9)./deltaTime'/1000;
%         bioTreeFrame=[];
%         nameList=dir(dirFluoFile);
%         for i=1:size(frameInfo,1)
%             bioTreeTimer=etime(frameInfo(i,1:6),bioTree{1}.bioTreeInitialTime)/60;
%             diffTime=abs(bioTree{1}.bioTreeTimer-bioTreeTimer);
%             [~,index]=find(diffTime==min(diffTime));
%             bioTreeFrame(i)=index(2);
%         end
%         [~,checkFrame]=find(diff(bioTreeFrame)==0);
%         bioTreeFrame(checkFrame+1)=0;
%         image=import_tiff_stack([dirFile,'\',fluoChannel{iChannel}(1:end-4),'Control\image',fluoChannel{iChannel}(1:end-4),'Control',num2str(1,'%05.f'),'.tif']);
%         for iPic=1:numel(bioTreeFrame)
%             if bioTreeFrame(iPic)==0
%                 continue
%             end
%             for iframe=1:bioTreeFrame(iPic)
%                 for iRoot=1:numel(bioTree{iframe}.root)
%                     traceInfo=bioTree{iframe}.root{iRoot}.traceInfo;
%                     bioTree{iframe}.root{iRoot}.linkTimer=bioTree{1}.bioTreeTimer(iframe:iframe+numel(traceInfo.pixelIdxList)-1)';
%                     if numel(traceInfo.pixelIdxList)+iframe>bioTreeFrame(iPic)
%                         tracePixel=traceInfo.pixelIdxList{bioTreeFrame(iPic)-iframe+1};
%                         switch iChannel
%                             case 8
%                                 bioTree{iframe}.root{iRoot}.RedControl(bioTreeFrame(iPic)-iframe+1)=mean(double(image(tracePixel)))*16/frameInfo(iPic,14)*frameInfo(iPic,15);
%                             case 9
%                                 bioTree{iframe}.root{iRoot}.BlueControl(bioTreeFrame(iPic)-iframe+1)=mean(double(image(tracePixel)))*16/frameInfo(iPic,14)*frameInfo(iPic,15);
%                             case 10
%                                 bioTree{iframe}.root{iRoot}.GreenControl(bioTreeFrame(iPic)-iframe+1)=mean(double(image(tracePixel)))*16/frameInfo(iPic,14)*frameInfo(iPic,15);
%                         end
%                     end
%                 end
%                 for iNode=1:numel(bioTree{iframe}.node)
%                     for iOut=1:numel(bioTree{iframe}.node{iNode}.Out)
%                         traceInfo=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo;
%                         bioTree{iframe}.node{iNode}.Out{iOut}.linkTimer=bioTree{1}.bioTreeTimer(iframe:iframe+numel(traceInfo.pixelIdxList)-1)';
%                         if numel(traceInfo.pixelIdxList)+iframe>bioTreeFrame(iPic)
%                             tracePixel=traceInfo.pixelIdxList{bioTreeFrame(iPic)-iframe+1};
%                             switch iChannel
%                                 case 8
%                                     bioTree{iframe}.node{iNode}.Out{iOut}.RedControl(bioTreeFrame(iPic)-iframe+1)=mean(double(image(tracePixel)))*16/frameInfo(iPic,14)*frameInfo(iPic,15);
%                                 case 9
%                                     bioTree{iframe}.node{iNode}.Out{iOut}.BlueControl(bioTreeFrame(iPic)-iframe+1)=mean(double(image(tracePixel)))*16/frameInfo(iPic,14)*frameInfo(iPic,15);
%                                 case 10
%                                     bioTree{iframe}.node{iNode}.Out{iOut}.GreenControl(bioTreeFrame(iPic)-iframe+1)=mean(double(image(tracePixel)))*16/frameInfo(iPic,14)*frameInfo(iPic,15);
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%         bioTree=bioTreeTraceConfirm(bioTree);    %   用前面的光强填充后面空白的部分
%     end
%     if strcmp(mip.Calib.lightControl,'dynamic') 
%         dirFluoFile=[dirFile,'\',fluoChannel{iChannel}];
%         try
%             load([dirFluoFile,'\frameInfo.mat'])
%         catch err
%             return
%         end
%         for i=1:size(frameInfo,1)-1
%             deltaTime(i)=abs(etime(frameInfo(i+1,1:6),frameInfo(i,1:6)));
%         end
%         deltaTime(end+1)=deltaTime(end);
%         frameInfo(:,15)=frameInfo(:,15).*frameInfo(:,9)./deltaTime'/1000;
%         bioTreeFrame=[];
%         nameList=dir(dirFluoFile);
%         for i=1:size(frameInfo,1)
%             bioTreeTimer=etime(frameInfo(i,1:6),bioTree{1}.bioTreeInitialTime)/60;
%             diffTime=abs(bioTree{1}.bioTreeTimer-bioTreeTimer);
%             [~,index]=find(diffTime==min(diffTime));
%             bioTreeFrame(i)=index(2);
%         end
%         [~,checkFrame]=find(diff(bioTreeFrame)==0);
%         bioTreeFrame(checkFrame+1)=0;
%         for iPic=1:numel(bioTreeFrame)
%             if bioTreeFrame(iPic)==0
%                 continue
%             end
%             image=import_tiff_stack([dirFile,'\',fluoChannel{iChannel}(1:end-4),'Control\image',fluoChannel{iChannel}(1:end-4),'Control',num2str(iPic,'%05.f'),'.tif']);
%             se=translate(strel(1),bestPositionAccumulation(bioTreeFrame(iPic),:));
%             image=imdilate(image,se);
%             for iframe=1:bioTreeFrame(iPic)
%                 for iRoot=1:numel(bioTree{iframe}.root)
%                     traceInfo=bioTree{iframe}.root{iRoot}.traceInfo;
%                     bioTree{iframe}.root{iRoot}.linkTimer=bioTree{1}.bioTreeTimer(iframe:iframe+numel(traceInfo.pixelIdxList)-1)';
%                     if numel(traceInfo.pixelIdxList)+iframe>bioTreeFrame(iPic)
%                         tracePixel=traceInfo.pixelIdxList{bioTreeFrame(iPic)-iframe+1};
%                         switch iChannel
%                             case 8
%                                 bioTree{iframe}.root{iRoot}.RedControl(bioTreeFrame(iPic)-iframe+1)=mean(double(image(tracePixel)))*16/frameInfo(iPic,14)*frameInfo(iPic,15);
%                             case 9
%                                 bioTree{iframe}.root{iRoot}.BlueControl(bioTreeFrame(iPic)-iframe+1)=mean(double(image(tracePixel)))*16/frameInfo(iPic,14)*frameInfo(iPic,15);
%                             case 10
%                                 bioTree{iframe}.root{iRoot}.GreenControl(bioTreeFrame(iPic)-iframe+1)=mean(double(image(tracePixel)))*16/frameInfo(iPic,14)*frameInfo(iPic,15);
%                         end
%                     end
%                 end
%                 for iNode=1:numel(bioTree{iframe}.node)
%                     for iOut=1:numel(bioTree{iframe}.node{iNode}.Out)
%                         traceInfo=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo;
%                         bioTree{iframe}.node{iNode}.Out{iOut}.linkTimer=bioTree{1}.bioTreeTimer(iframe:iframe+numel(traceInfo.pixelIdxList)-1)';
%                         if numel(traceInfo.pixelIdxList)+iframe>bioTreeFrame(iPic)
%                             tracePixel=traceInfo.pixelIdxList{bioTreeFrame(iPic)-iframe+1};
%                             switch iChannel
%                                 case 8
%                                     bioTree{iframe}.node{iNode}.Out{iOut}.RedControl(bioTreeFrame(iPic)-iframe+1)=mean(double(image(tracePixel)))*16/frameInfo(iPic,14)*frameInfo(iPic,15);
%                                 case 9
%                                     bioTree{iframe}.node{iNode}.Out{iOut}.BlueControl(bioTreeFrame(iPic)-iframe+1)=mean(double(image(tracePixel)))*16/frameInfo(iPic,14)*frameInfo(iPic,15);
%                                 case 10
%                                     bioTree{iframe}.node{iNode}.Out{iOut}.GreenControl(bioTreeFrame(iPic)-iframe+1)=mean(double(image(tracePixel)))*16/frameInfo(iPic,14)*frameInfo(iPic,15);
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%         bioTree=bioTreeTraceConfirm(bioTree);    %   用前面的光强填充后面空白的部分
%     end
% end
end
function bioTree=bioTreeAllChannel(bioTree)
for iframe=1:numel(bioTree)
    for iRoot=1:numel(bioTree{iframe}.root)
        bioTree{iframe}.root{iRoot}.CyOFP=ones(numel(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList),1)*nan;
        bioTree{iframe}.root{iRoot}.GFP=ones(numel(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList),1)*nan;
        bioTree{iframe}.root{iRoot}.mScalet=ones(numel(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList),1)*nan;
        bioTree{iframe}.root{iRoot}.RFP=ones(numel(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList),1)*nan;
        bioTree{iframe}.root{iRoot}.CyPet=ones(numel(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList),1)*nan;
        bioTree{iframe}.root{iRoot}.Venus=ones(numel(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList),1)*nan;
        bioTree{iframe}.root{iRoot}.mAmetrine=ones(numel(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList),1)*nan;
        bioTree{iframe}.root{iRoot}.RedControl=ones(numel(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList),1)*0;
        bioTree{iframe}.root{iRoot}.BlueControl=ones(numel(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList),1)*0;
        bioTree{iframe}.root{iRoot}.GreenControl=ones(numel(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList),1)*0;
    end
    for iNode=1:numel(bioTree{iframe}.node)
        for iOut=1:numel(bioTree{iframe}.node{iNode}.Out)
            bioTree{iframe}.node{iNode}.Out{iOut}.CyOFP=ones(numel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList),1)*nan;
            bioTree{iframe}.node{iNode}.Out{iOut}.GFP=ones(numel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList),1)*nan;
            bioTree{iframe}.node{iNode}.Out{iOut}.mScalet=ones(numel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList),1)*nan;
            bioTree{iframe}.node{iNode}.Out{iOut}.RFP=ones(numel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList),1)*nan;
            bioTree{iframe}.node{iNode}.Out{iOut}.CyPet=ones(numel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList),1)*nan;
            bioTree{iframe}.node{iNode}.Out{iOut}.Venus=ones(numel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList),1)*nan;
            bioTree{iframe}.node{iNode}.Out{iOut}.mAmetrine=ones(numel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList),1)*nan;
            bioTree{iframe}.node{iNode}.Out{iOut}.RedControl=ones(numel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList),1)*0;
            bioTree{iframe}.node{iNode}.Out{iOut}.BlueControl=ones(numel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList),1)*0;
            bioTree{iframe}.node{iNode}.Out{iOut}.GreenControl=ones(numel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList),1)*0;
        end
    end
end
end
function bioTree=bioTreeTraceConfirm(bioTree)
for iframe=1:numel(bioTree)
    for iRoot=1:numel(bioTree{iframe}.root)
        bioTree{iframe}.root{iRoot}.RedControl=selfFillData(bioTree{iframe}.root{iRoot}.RedControl);
        bioTree{iframe}.root{iRoot}.BlueControl=selfFillData(bioTree{iframe}.root{iRoot}.BlueControl);
        bioTree{iframe}.root{iRoot}.GreenControl=selfFillData(bioTree{iframe}.root{iRoot}.GreenControl);
    end
    for iNode=1:numel(bioTree{iframe}.node)
        for iIn=1:numel(bioTree{iframe}.node{iNode}.In)
            if bioTree{iframe}.node{iNode}.In{iIn}.isNode==0
                rootInfo=bioTree{iframe}.node{iNode}.In{iIn}.rootInfo;
                RedlightSource=bioTree{rootInfo(1)}.root{rootInfo(2)}.RedControl(end);
                BluelightSource=bioTree{rootInfo(1)}.root{rootInfo(2)}.BlueControl(end);
                GreenlightSource=bioTree{rootInfo(1)}.root{rootInfo(2)}.GreenControl(end);
            end
            if bioTree{iframe}.node{iNode}.In{iIn}.isNode==1
                nodeInfo=bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo;
                RedlightSource=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.RedControl(end);
                BluelightSource=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.BlueControl(end);
                GreenlightSource=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.GreenControl(end);
            end
        end
        for iOut=1:numel(bioTree{iframe}.node{iNode}.Out)
            if bioTree{iframe}.node{iNode}.Out{iOut}.RedControl(1)==0
                bioTree{iframe}.node{iNode}.Out{iOut}.RedControl=max(RedlightSource);
            end
            if bioTree{iframe}.node{iNode}.Out{iOut}.BlueControl(1)==0
                bioTree{iframe}.node{iNode}.Out{iOut}.BlueControl=max(BluelightSource);
            end
            if bioTree{iframe}.node{iNode}.Out{iOut}.GreenControl(1)==0
                bioTree{iframe}.node{iNode}.Out{iOut}.GreenControl=max(GreenlightSource);
            end
            bioTree{iframe}.node{iNode}.Out{iOut}.RedControl=selfFillData(bioTree{iframe}.node{iNode}.Out{iOut}.RedControl);
            bioTree{iframe}.node{iNode}.Out{iOut}.BlueControl=selfFillData(bioTree{iframe}.node{iNode}.Out{iOut}.BlueControl);
            bioTree{iframe}.node{iNode}.Out{iOut}.GreenControl=selfFillData(bioTree{iframe}.node{iNode}.Out{iOut}.GreenControl);
        end
    end
end
end
function dataLine=selfFillData(dataLine)
if ~all((dataLine)==0)
    index=find(dataLine~=0);
    if index(end)~=numel(dataLine)
        index(end+1)=numel(dataLine);
    end
    for i=1:numel(index)-1
        dataLine(index(i):index(i+1)-1)=dataLine(index(i));
    end
end
end