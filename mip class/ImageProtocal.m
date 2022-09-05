classdef ImageProtocal
    
    properties
        ImagePara = ImageParameter;
        DeviceS = DeviceState;
        Calib = Calibration;
        xyzC = xyzControl;
        currentImage = currentImage;
        imageProcessing = imageProcessing;
        hardware = hardwareSetup;
        fieldTag = fieldTag;
        niDevices;
        nidaq; %nidaq
        niSignal;
        imageShowing = imageShowing;
        gui %micromanager gui
        mmc_IX83 %Olympus IX83 control
        mmc_ASI %XY stage and piezo Z control
        mmc_Camera1 %Prime BSI camera 1
        experimentIdx;
        myDir
        experimentType;
        protocalName;
    end
    
    methods
        function obj = initialization(obj)
            %when run the function,code must be 'obj = obj.initialization();'. It is necessary to have a returned value.If not, mmc will be losed and cannot unload the devices.
            %             addpath('xxxxx')                          %
            obj.myDir = uigetdir();
            %             mmc=loadMMconfigFile('C:\Program Files\Micro-Manager-1.4.23\allDevices2.cfg');
            try
                disp('Please waiting for 30s. Connect to all Devices');
                obj.niDevices = daq.getDevices;
                obj.nidaq = daq.createSession('ni') ;
                addDigitalChannel(obj.nidaq,'PXI1Slot5','port0/line0','OutputOnly');
                addDigitalChannel(obj.nidaq,'PXI1Slot5','port0/line1','OutputOnly');
                addAnalogOutputChannel(obj.nidaq,'PXI1Slot5','ao0','Voltage');
                obj.nidaq.Rate = 100000; %Will run at 100000 scans/second.
                lh = addlistener(obj.nidaq,'DataAvailable',@plotData);
                obj.DeviceS.PiezoZVoltage = 5;
                outputSingleScan(obj.nidaq,[0 0 5]);%give an initial voltage, camera 0 V,led 0 V,piezoZ 5 V
                
                %                 obj.mmc = loadMMconfigFile('C:\Program Files\Micro-Manager-2.0gamma\20191112_XAG_BSI1_ASIXY.cfg');
                obj.gui = StartMMStudio();
                
%                 obj.mmc = loadMMconfigFile('C:\Program Files\Micro-Manager-2.0gamma\20191121_XAG_BSI1_ASIXYZ_IX83.cfg');
                obj.mmc_IX83 = loadMMconfigFile('C:\Program Files\Micro-Manager-2.0gamma\20191130_XAG_IX83.cfg');
                obj.mmc_ASI = loadMMconfigFile('C:\Program Files\Micro-Manager-2.0gamma\20191130_XAG_ASI.cfg');
                obj.mmc_Camera1 = obj.gui.core;
                disp('Import Fiji..');
                addpath('D:\Fiji\Fiji.app\scripts');
                Miji;
                import ij.IJ.*;
                disp('Setup hareware...');
                obj.mmc_Camera1.setConfig('channel','BF');
                %                 obj.mmc.setConfig('Olympus','IX83');
                %                 obj.mmc.setConfig('ZDC','Off');
                disp('Welcome to a new experiment');
                obj = hardwareSet(obj);
                obj.currentImage = obj.currentImage.initializeCurrentImage();
                obj.imageShowing = obj.imageShowing.initializeFigureSetting();
                obj.Calib = obj.Calib.initializeCalibration();
                obj.experimentIdx = input('Please input experimental idx:');
                obj.experimentType = input('Please input experimental type:1 flow cell, 2 static cell, 3 plate');
                olympus = input('Which light path of olympus will you select:''1x'' or ''1.6x'';');
                lightControl = input('Please input mode of light control:''no'',''static'',''dynamic'';');
                obj.Calib.Olympus = olympus;
                obj.Calib.lightControl = lightControl;
                objective = input('Which objective of olympus will you select:[100];');
                obj.Calib.Objective = objective;
                switch olympus
                    case '1x'
                        obj.Calib.scale = 0.065*100/objective;
                    case '1.6x'
                        obj.Calib.scale = 0.065/1.6*100/objective;
                end
                %                 if obj.Calib.Objective==20
                %                     obj.mmc.setProperty('CRISP','Calibration Gain','80');
                %                     obj.mmc.setProperty('CRISP','GainMultiplier','10');
                %                     obj.mmc.setProperty('CRISP','Objective NA','0.86');
                %                 end
            catch ME
                obj.mmc_IX83.unloadAllDevices();
                obj.mmc_ASI.unloadAllDevices();
                obj.mmc_Camera1.unloadAllDevices();
%                 obj.mmc.unloadAllDevices();
                obj.nidaq.removeChannel(3);
                obj.nidaq.removeChannel(2);
                obj.nidaq.removeChannel(1);
                error('Initialization failed');
            end
        end
        
        function obj = hardwareSet(obj)
            
            obj.mmc_Camera1.setCircularBufferMemoryFootprint(5000);%Reserve memory for the circular buffer.
            
            obj.hardware.PhysicalCamera1Trigger = char(obj.mmc_Camera1.getProperty('Physical Camera 1','TriggerMode'));
            obj.hardware.PhysicalCamera1Exposure = char(obj.mmc_Camera1.getProperty('Physical Camera 1','Exposure'));
            obj.hardware.PhysicalCamera1ReadoutRate = char(obj.mmc_Camera1.getProperty('Physical Camera 1','ReadoutRate'));
            %             obj.hardware.CoolLedIntensityA = char(obj.mmc.getProperty('pE4000','IntensityA'));
            %             obj.hardware.CoolLedIntensityB = char(obj.mmc.getProperty('pE4000','IntensityB'));
            %             obj.hardware.CoolLedIntensityC = char(obj.mmc.getProperty('pE4000','IntensityC'));
            %             obj.hardware.CoolLedIntensityD = char(obj.mmc.getProperty('pE4000','IntensityD'));
            %             obj.hardware.CoolLedSelectionA = char(obj.mmc.getProperty('pE4000','SelectionA'));
            %             obj.hardware.CoolLedSelectionB = char(obj.mmc.getProperty('pE4000','SelectionB'));
            %             obj.hardware.CoolLedSelectionC = char(obj.mmc.getProperty('pE4000','SelectionC'));
            %             obj.hardware.CoolLedSelectionD = char(obj.mmc.getProperty('pE4000','SelectionD'));
        end
        
        function obj = getCurrentDeviceState(obj)
            obj.DeviceS.PhysicalCamera1Trigger = char(obj.mmc_Camera1.getProperty('Physical Camera 1','TriggerMode'));
            obj.DeviceS.PhysicalCamera1Exposure = char(obj.mmc_Camera1.getProperty('Physical Camera 1','Exposure'));
            obj.DeviceS.PhysicalCamera1ReadoutRate = char(obj.mmc_Camera1.getProperty('Physical Camera 1','ReadoutRate'));
            obj.DeviceS.CameraIdx = 'Physical Camera 1';
            %             obj.DeviceS.CoolLedIntensityA = char(obj.mmc.getProperty('pE4000','IntensityA'));
            %             obj.DeviceS.CoolLedIntensityB = char(obj.mmc.getProperty('pE4000','IntensityB'));
            %             obj.DeviceS.CoolLedIntensityC = char(obj.mmc.getProperty('pE4000','IntensityC'));
            %             obj.DeviceS.CoolLedIntensityD = char(obj.mmc.getProperty('pE4000','IntensityD'));
            %             obj.DeviceS.CoolLedSelectionA = char(obj.mmc.getProperty('pE4000','SelectionA'));
            %             obj.DeviceS.CoolLedSelectionB = char(obj.mmc.getProperty('pE4000','SelectionB'));
            %             obj.DeviceS.CoolLedSelectionC = char(obj.mmc.getProperty('pE4000','SelectionC'));
            %             obj.DeviceS.CoolLedSelectionD = char(obj.mmc.getProperty('pE4000','SelectionD'));
            obj.DeviceS.PiezoZVoltage = 5;
            
        end
        
        function obj = setCameraIdx(obj,CameraIdx) %set only first time
            if ~strcmp(obj.DeviceS.CameraIdx,CameraIdx)
                obj.mmc_Camera1.setCameraDevice(CameraIdx);
                obj.DeviceS.CameraIdx = char(obj.mmc_Camera1.getCameraDevice());
            end
            switch CameraIdx
                case 'Physical Camera 1'
                    if ~strcmp(char(obj.DeviceS.PhysicalCamera1Trigger),'Internal Trigger')
                        obj.mmc_Camera1.setProperty('Physical Camera 1','TriggerMode','Internal Trigger'); %Trigger List:Internal Trigger,Edge Trigger
                        obj.DeviceS.PhysicalCamera1Trigger = 'Internal Trigger';
                    end
                    if ~strcmp(char(obj.DeviceS.PhysicalCamera2Trigger),'Internal Trigger')
                        obj.mmc_Camera1.setProperty('Physical Camera 2','TriggerMode','Internal Trigger');
                        obj.DeviceS.PhysicalCamera2Trigger = 'Internal Trigger';
                    end
                case 'Physical Camera 2'
                    if ~strcmp(char(obj.DeviceS.PhysicalCamera1Trigger),'Internal Trigger')
                        obj.mmc_Camera1.setProperty('Physical Camera 1','TriggerMode','Internal Trigger');
                        obj.DeviceS.PhysicalCamera1Trigger = 'Internal Trigger';
                    end
                    if ~strcmp(char(obj.DeviceS.PhysicalCamera2Trigger),'Internal Trigger')
                        obj.mmc_Camera1.setProperty('Physical Camera 2','TriggerMode','Internal Trigger');
                        obj.DeviceS.PhysicalCamera2Trigger = 'Internal Trigger';
                    end
                case 'Multi Camera'
                    if ~strcmp(char(obj.DeviceS.PhysicalCamera1Trigger),'Internal Trigger')
                        obj.mmc_Camera1.setProperty('Physical Camera 1','TriggerMode','Internal Trigger');
                        obj.DeviceS.PhysicalCamera1Trigger = 'Internal Trigger';
                    end
                    if ~strcmp(char(obj.DeviceS.PhysicalCamera2Trigger),'Edge Trigger')
                        obj.mmc_Camera1.setProperty('Physical Camera 2','TriggerMode','Edge Trigger');
                        obj.DeviceS.PhysicalCamera2Trigger = 'Edge Trigger';
                    end
            end
        end
        
        
        function obj=setChannelDeviceState(obj,channel)
            switch channel
                case 'corr_BF_hardware'
                    
                    mip.mmc_Camera1.setProperty('Physical Camera 1','TriggerMode','Edge Trigger');%select camera trigger mode to edge trigger;

                    mip.mmc_Camera1.setProperty('Physical Camera 1','ReadoutRate','100MHz 16bit');
                case 'BF'
                    mip.mmc_Camera1.setProperty('Physical Camera 1','TriggerMode','Internal Trigger');
                    %                     obj.mmc.setProperty('Arduino-Switch','Label','1.ch1');%select BF LED;
                    %                     obj.mmc.setProperty('Arduino-Shutter','OnOff','0');
                    %                     if ~strcmp(obj.DeviceS.FilterCube,'Position-3')
                    %                         obj.mmc.setProperty('FilterCube','Label','Position-3');
                    %                         obj.mmc.waitForDevice('FilterCube');
                    %                         obj.DeviceS.FilterCube = char(obj.mmc.getProperty('FilterCube','Label'));
                    %                     end
                    %                     if ~strcmp(obj.DeviceS.WheelA,'Filter-1')
                    %                         obj.mmc.setProperty('Wheel-A','Label','Filter-1');
                    %                         obj.mmc.waitForDevice('Wheel-A');
                    %                         obj.DeviceS.WheelA =char(obj.mmc.getProperty('Wheel-A','Label'));
                    %                     end
                    %                     if ~strcmp(obj.DeviceS.WheelB,'Filter-1')
                    %                         obj.mmc.setProperty('Wheel-B','Label','Filter-1');
                    %                         obj.mmc.waitForDevice('Wheel-B');
                    %                         obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    %                     end
                    %                     if ~strcmp(obj.DeviceS.Shutter,'Arduino-Shutter')
                    %                         obj.mmc.setShutterDevice('Arduino-Shutter');
                    %                         obj.mmc.waitForDevice('Core');
                    %                         obj.DeviceS.Shutter =char(obj.mmc.getShutterDevice());
                    %                     end
%                     if obj.mmc_Camera1.getExposure ~=obj.ImagePara.BFexposeTime
%                         obj.mmc_Camera1.setExposure(obj.ImagePara.BFexposeTime);
%                     end
            end
        end
        
        function [obj,imageStack]=sequenceAcquisitionImage(obj,nrFrames,cameraVoltage,ledVoltage,piezoVoltage)
            width = obj.currentImage.imageSize(1);
            height = obj.currentImage.imageSize(2);
            obj.mmc_Camera1.initializeCircularBuffer();%Initialize circular buffer based on the current camera settings.
            % gui.initializeAcquisition(acqName, width, height, bytesPerPixel, bitDepth);
            %             obj.mmc.prepareSequenceAcquisition('Physical Camera 1');
            % frame = 0;
            exposureMs = obj.mmc_Camera1.getExposure();
            camera_bit = char(obj.mmc_Camera1.getProperty('Physical Camera 1','ReadoutRate'));
            switch camera_bit
                case '200MHz 11bit'
                    obj.mmc_Camera1.setExposure(obj.ImagePara.BFexposeTime + 16);% Sets the exposure setting of the current camera in milliseconds.
                case '100MHz 16bit'
                    obj.mmc_Camera1.setExposure(obj.ImagePara.BFexposeTime + 24);% Sets the exposure setting of the current camera in milliseconds.
            end
            %             obj.mmc.setExposure(fix(numel(cameraVoltage)/obj.nidaq.Rate*1000/2-1));% Sets the exposure setting of the current camera in milliseconds.
            %             tic;
            obj.mmc_Camera1.prepareSequenceAcquisition('Physical Camera 1');
            %             mmcSlaveCam.prepareSequenceAcquisition('Physical Camera 2');
            % mmc.prepareSequenceAcquisition('Multi Camera');%时间已过 3.631730 秒。
            % tic;
            disp('start image acquisition');
            obj.mmc_Camera1.startSequenceAcquisition(nrFrames, 0, true);%(int numImages, double intervalMs, boolean stopOnOverflow)
            %             mmcSlaveCam.startSequenceAcquisition(nrFrames, 0, true);%(int numImages, double intervalMs, boolean stopOnOverflow)
            % toc;
            %             disp(obj.mmc.getRemainingImageCount());
            if min(piezoVoltage) >= 0 && max(piezoVoltage) <= 10
                %     tic;queueOutputData(s,[step1(1:1000000),step3(1:1000000),step2(1:1000000)]);
                %                 tic;
                queueOutputData(obj.nidaq,[cameraVoltage,ledVoltage,piezoVoltage]);
                obj.nidaq.startBackground();
                %                 toc;% 0.101141s
            else
                obj.mmc_Camera1.stopSequenceAcquisition();%Stops streaming camera sequence acquisition.
                disp('External voltage of PizeoZ must be 0~10 V!');
                return;
            end
            imageStack = zeros(width*height*nrFrames,1,'int16');
            % imageStack1=[];
            %             imageStack2 = zeros(width*height*nrFrames,1,'int16');
            %             tic
            % spmd
            %     if labindex ==2
            frame = 0;
            
            while ((obj.mmc_Camera1.getRemainingImageCount() > 0) || (obj.mmc_Camera1.isSequenceRunning(obj.mmc_Camera1.getCameraDevice())))
                if (obj.mmc_Camera1.getRemainingImageCount() > 0)
                    %                     disp(obj.mmc.getRemainingImageCount());
                    imageStack(frame*width*height+1:(frame+1)*width*height) = obj.mmc_Camera1.popNextImage(); %tic;toc;62 frame 1.34s
                    frame = frame + 1;
                    %    else
                    %       pause(1);
                end
            end
            disp(['get all images ' num2str(frame)]);
            % toc
            imageStack = reshape(imageStack,[width,height,nrFrames]);
            %             toc
            %             imageStack2=reshape(imageStack2,[width,height,nrFrames]);
            %             toc
        end
        
        function obj= getSnapImage(obj,channelIdx)
%             camera_bit = char(obj.mmc.getProperty('Physical Camera 1','ReadoutRate'));
%             switch camera_bit
%                 case '200MHz 11bit'
%                     obj.mmc.setExposure(obj.ImagePara.BFexposeTime+16);% Sets the exposure setting of the current camera in milliseconds.
%                 case '100MHz 16bit'
%                     obj.mmc.setExposure(obj.ImagePara.BFexposeTime+24);% Sets the exposure setting of the current camera in milliseconds.
%             end
            switch obj.DeviceS.CameraIdx
                case 'Physical Camera 1'
                    switch channelIdx
                        case 'BF'
                            %                     obj.mmc.setProperty('Arduino-Shutter','OnOff','1');%BF LED on
                            obj.mmc_Camera1.snapImage(); %after set shutter to Arduino-Switch,run mmc.snapImage will tigger BF LED on
%                             outputSingleScan(obj.nidaq,[1 0 5]);
%                             switch camera_bit
%                                 case '200MHz 11bit'
%                                     pause(16/1000);% Sets the exposure setting of the current camera in milliseconds.
%                                 case '100MHz 16bit'
%                                     pause(24/1000);% Sets the exposure setting of the current camera in milliseconds.
%                             end                           
%                             outputSingleScan(obj.nidaq,[0 1 5]);
%                             pause(mip.ImagePara.BFexposeTime/1000);
%                             outputSingleScan(obj.nidaq,[0 0 5]);
                            %                             obj.mmc.setProperty('Arduino-Shutter','OnOff','0');%BF LED off, setting Arduino-Shutter off will faster turnoff the BF LED
                        case 'Mosaic'
                            obj.mmc_Camera1.displaySLMImage(obj.DeviceS.Mosaic);
                            obj.mmc_Camera1.snapImage();
                        otherwise
                            disp('None channelIdx');
                    end
                    
                    image1 = obj.mmc_Camera1.getImage(0);  % Physical Camera 1 get image
                    %                     image2 = obj.mmc.getImage(1);  % Physical Camera 2 get image
                    width = obj.mmc_Camera1.getImageWidth();
                    height = obj.mmc_Camera1.getImageHeight();
                    if obj.mmc_Camera1.getBytesPerPixel == 2
                        pixelType = 'uint16';
                    else
                        pixelType = 'uint8';
                    end
                    image1 = typecast(image1, pixelType); % pixels must be interpreted as unsigned integers
                    image1 = reshape(image1, [width, height]); % image should be interpreted as a 2D array
                    image1 = transpose(image1);
                    %two camera adjustion must
                    image1 = flipud(image1);%image of Physical Camera 1 is now up and down symmetry.If images are left and right symmetry,use function 'fliplr'
                    
                    %                     image2 = typecast(image2, pixelType); % pixels must be interpreted as unsigned integers
                    %                     image2 = reshape(image2, [width, height]); % image should be interpreted as a 2D array
                    %                     image2 = transpose(image2);
                    %
                    obj.currentImage.imageOri(:,:,1) = image1;
                    %                     obj.currentImage.imageOri(:,:,2)=image2;
                case 'Physical Camera 2'
                    
                    
                case 'Multi Camera'
                    switch channelIdx
                        case 'BF'
                            %                     obj.mmc.setProperty('Arduino-Shutter','OnOff','1');%BF LED on
                            obj.mmc_Camera1.snapImage(); %after set shutter to Arduino-Switch,run mmc.snapImage will tigger BF LED on
                            obj.mmc_Camera1.setProperty('Arduino-Shutter','OnOff','0');%BF LED off, setting Arduino-Shutter off will faster turnoff the BF LED
                        case 'Mosaic'
                            obj.mmc_Camera1.displaySLMImage(obj.DeviceS.Mosaic);
                            obj.mmc_Camera1.snapImage();
                        otherwise
                            disp('None channelIdx');
                    end
                    
                    image1 = obj.mmc_Camera1.getImage(0);  % Physical Camera 1 get image
                    image2 = obj.mmc_Camera1.getImage(1);  % Physical Camera 2 get image
                    width = obj.mmc_Camera1.getImageWidth();
                    height = obj.mmc_Camera1.getImageHeight();
                    if obj.mmc_Camera1.getBytesPerPixel == 2
                        pixelType = 'uint16';
                    else
                        pixelType = 'uint8';
                    end
                    image1 = typecast(image1, pixelType); % pixels must be interpreted as unsigned integers
                    image1 = reshape(image1, [width, height]); % image should be interpreted as a 2D array
                    image1 = transpose(image1);
                    %two camera adjustion must
                    image1=flipud(image1);%image of Physical Camera 1 is now up and down symmetry.If images are left and right symmetry,use function 'fliplr'
                    
                    image2 = typecast(image2, pixelType); % pixels must be interpreted as unsigned integers
                    image2 = reshape(image2, [width, height]); % image should be interpreted as a 2D array
                    image2 = transpose(image2);
                    
                    obj.currentImage.imageOri(:,:,1)=image1;
                    obj.currentImage.imageOri(:,:,2)=image2;
                    
            end
            
        end
        
        
        function obj=imageshow(obj,varagin)
            imageNum=numel(varagin);
            for i=1:imageNum
                switch varagin{i}
                    case 'imageBFOri'
                        figure;
                        if ~isempty(obj.currentImage.imageBFOri)
                            imshowpair(imadjust(obj.currentImage.imageBFOri(:,:,1)),imadjust(obj.currentImage.imageBFOri(:,:,2)),'montage');
                        else
                            disp('imageBFOri is empty');
                        end
                    case 'imageBF'
                        figure;
                        if ~isemtpy(obj.currentImage.imageBF)
                            imshowpair(imadjust(obj.currentImage.imageBF(:,:,1)),imadjust(obj.currentImage.imageBF(:,:,2)),'montage');
                        else
                            disp('imageBFOri is empty');
                        end
                    case 'imageFluoOri'
                        figure;
                        if ~isempty(obj.currentImage.imageFluoOri)
                            imshowpair(imadjust(obj.currentImage.imageFluoOri(:,:,1)),imadjust(obj.currentImage.imageFluoOri(:,:,2)),'montage');
                        else
                            disp('imageFluoOri is empty');
                        end
                    case 'imageFluo'
                        figure;
                        if ~isempty(obj.currentImage.imageFluo)
                            imshowpair(imadjust(obj.currentImage.imageFluo(:,:,1)),imadjust(obj.currentImage.imageFluo(:,:,1)),'montage');
                        else
                            disp('imageFluo is empty');
                        end
                end
                
                
                
            end
        end
        
%         function obj = 
%         
%         end
        
        function obj = setZPostion(obj,offsetZ)
            error=0.15;
            currentPosition=obj.mmc.getPosition();
            %             disp(num2str(obj.mmc.getPosition()));
            obj.mmc.setPosition(currentPosition-5);
            obj.mmc.waitForDevice('ZStage');
            %             disp(num2str(obj.mmc.getPosition()));
            zPosition=currentPosition+offsetZ;
            %             while abs(zPosition-currentPosition)>error
            %                 offsetZ=zPosition-currentPosition;
            obj.mmc.setPosition(zPosition)
            obj.mmc.waitForDevice('ZStage');
            %             disp(num2str(obj.mmc.getPosition()));
            %             currentPosition=obj.mmc.getPosition();
            %             end
        end
        
        
        
        function obj = getCalibrationData(obj,lightSelect)
            calibfieldNum=size(obj.Calib.calibrationFeild,1);
            obj.Calib.calibrationFeild(calibfieldNum+1,1)=obj.mmc_ASI.getXPosition();
            obj.Calib.calibrationFeild(calibfieldNum+1,2)=obj.mmc_ASI.getYPosition();
            obj=obj.getCurrentDeviceState();
            if ~isempty(lightSelect)
                for i=1:numel(lightSelect)
                    switch lightSelect{i}
                        case 'BlueControl'
                            %                             obj = obj.moveFocalPlane('BlueControl',1);
                            obj=obj.setChannelDeviceState('BlueControl');
                            [obj.Calib.transBlueControl,~,~]=imageTransformMosaic2Zyla(obj.mmc,100,{'Blue'},{'10'});
                        case 'GreenControl'
                            obj = obj.moveFocalPlane('GreenControl',0);
                            obj=obj.setChannelDeviceState('GreenControl');
                            [obj.Calib.transGreenControl,~,~]=imageTransformMosaic2Zyla(obj.mmc,100,{'Teal'},{'80'});
                        case 'RedControl'
                            obj = obj.moveFocalPlane('RedControl',0);
                            obj=obj.setChannelDeviceState('RedControl');
                            [obj.Calib.transRedControl,~,~]=imageTransformMosaic2Zyla(obj.mmc,50,{'Red'},{'10'});
                    end
                end
            end
        end
        
        function obj=addXYPosition(obj)
            
            if isempty(obj.fieldTag.name)
                obj.fieldTag.name{1} = 'field0001';
                obj.fieldTag.xyPosition(1,1) = obj.mmc_ASI.getXPosition();
                obj.fieldTag.xyPosition(1,2) = obj.mmc_ASI.getYPosition();
                obj.fieldTag.tag{1} = 'none';
                obj.fieldTag.tagValue(1,:) = 0;
            else
                currentfieldTag = size(obj.fieldTag.name,2);
                obj.fieldTag.name{currentfieldTag+1} = strcat('field',num2str(currentfieldTag+1,'%.4d'));
                obj.fieldTag.xyPosition(currentfieldTag+1,1) = obj.mmc_ASI.getXPosition();
                obj.fieldTag.xyPosition(currentfieldTag+1,2) = obj.mmc_ASI.getYPosition();
                obj.fieldTag.tag{currentfieldTag+1} = 'none';
                obj.fieldTag.tagValue(currentfieldTag+1,:) = 0;
            end
        end
        
        function obj=addMontagefield(obj,xnumber,ynumber,gap)
            %             if ~isempty(obj.fieldTag.name)
            %                 obj.fieldTag.name={};
            %                 obj.fieldTag.xyPosition=[];
            %                 obj.fieldTag.tag={};
            %                 obj.fieldTag.tagValue=[];
            %             end
            currentfieldTag1 = size(obj.fieldTag.name,2);
            imageSize = [2048,2048];
            currentXYPosition(1) = obj.mmc_ASI.getXPosition();
            currentXYPosition(2) = obj.mmc_ASI.getYPosition();
            for x=1:xnumber
                for y=1:ynumber                    
                    currentfieldTag = (x-1)*ynumber+y+currentfieldTag1;
                    obj.fieldTag.name{currentfieldTag} = strcat('field',num2str(currentfieldTag,'%04d'));
                    obj.fieldTag.xyPosition(currentfieldTag,1:2) = [currentXYPosition(1)+(1+gap)*(x-1)*imageSize(1)*obj.Calib.scale,currentXYPosition(2)+(1+gap)*(y-1)*imageSize(2)*obj.Calib.scale];
                    obj.fieldTag.tag{currentfieldTag} = 'none';
                    obj.fieldTag.tagValue(currentfieldTag,:) = 0;
                end
            end
            
            
        end
        
        function obj=deleteXYPosition(obj,varargin)
            if isempty(varargin)
                obj.fieldTag.name={};
                obj.fieldTag.xyPosition=[];
                obj.fieldTag.tag={};
                obj.fieldTag.tagValue=[];
            else
                delectXYfield=varargin{1};
                delectXYfield=sort(delectXYfield,'descend');
                for ifield=1:numel(delectXYfield)
                    obj.fieldTag.name(delectXYfield(ifield))=[];
                    obj.fieldTag.xyPosition(delectXYfield(ifield),:)=[];
                    obj.fieldTag.tag(delectXYfield(ifield))=[];
                    obj.fieldTag.tagValue(delectXYfield(ifield))=[];
                    if delectXYfield(ifield)<size(obj.fieldTag.name,1)+1
                        for i=delectXYfield(ifield):size(obj.fieldTag.name,2)
                            obj.fieldTag.name{i}=strcat('field',num2str(i,'%.4d'));
                        end
                    end
                    
                end
            end
        end
        
        function obj=showfieldPosition(obj)
            imageSize=[2048,2048];
            rectangleW=obj.Calib.scale*imageSize(1);
            rectangleH=obj.Calib.scale*imageSize(2);
            figure;
            hold on;
            axis equal;
            if ~isempty(obj.Calib.calibrationFeild)
                for i=1:size(obj.Calib.calibrationFeild,1)
                    label='Cal';
                    plot(obj.Calib.calibrationFeild(i,1),obj.Calib.calibrationFeild(i,2));
                    text(obj.Calib.calibrationFeild(i,1),obj.Calib.calibrationFeild(i,2),label,'Color',[1,0,0]);
                    rectangle('Position',[obj.Calib.calibrationFeild(i,1)-rectangleW/2,obj.Calib.calibrationFeild(i,2)-rectangleH/2,rectangleW,rectangleH], 'LineWidth',2,...
                        'LineStyle','--',...
                        'EdgeColor',[1 0 0]);
                end
            end
            if ~isempty(obj.fieldTag.name)
                for i=1:size(obj.fieldTag.xyPosition,1)
                    label=num2str(i);
                    plot(obj.fieldTag.xyPosition(i,1),obj.fieldTag.xyPosition(i,2));
                    text(obj.fieldTag.xyPosition(i,1),obj.fieldTag.xyPosition(i,2),label,'Color',[1,0,0]);
                    rectangle('Position',[obj.fieldTag.xyPosition(i,1)-rectangleW/2,obj.fieldTag.xyPosition(i,2)-rectangleH/2,rectangleW,rectangleH], 'LineWidth',2,...
                        'LineStyle','--',...
                        'EdgeColor',[0 0.498039215803146 0]);
                end
            end
            
        end
        
        function obj=move2field(obj,fieldNum)
            obj.mmc_ASI.setXYPosition(obj.fieldTag.xyPosition(fieldNum,1),obj.fieldTag.xyPosition(fieldNum,2));
            %             obj.mmc.waitForDevice('XYStage');
            obj.xyzC.currentFeild=obj.fieldTag.name{fieldNum};
            obj.xyzC.XYpositions(1)=obj.mmc_ASI.getXPosition();
            obj.xyzC.XYpositions(2)=obj.mmc_ASI.getYPosition();
        end
        
        function obj=saveImageFile(obj,tag,imageProcessing)
            
            currentfield=obj.xyzC.currentFeild;
            currentfieldFile=strcat(obj.myDir,'\',currentfield);
            if ~exist(currentfieldFile,'file')
                mkdir(currentfieldFile);
            end
            feilTag=strcat(currentfieldFile,'\',tag);
            if ~exist(feilTag,'file')
                mkdir(feilTag);
            end
            
            if exist(strcat(feilTag,'\','frameInfo','.mat'),'file')
                load(strcat(feilTag,'\','frameInfo','.mat'));
                frameNum=size(frameInfo,1)+1;
            else
                frameNum=1;
                
            end
            
            frameInfo(frameNum,1:6)=clock;
            frameInfo(frameNum,7)=obj.mmc_ASI.getXPosition();
            frameInfo(frameNum,8)=obj.mmc_ASI.getYPosition();
            if imageProcessing==1
                frameInfo(frameNum,11)=obj.imageProcessing.cellNum;
                frameInfo(frameNum,12)=obj.imageProcessing.cellTotalArea;
            end
            save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
            %             timeInfo=eval(strcat(tag,'timeInfo'));
            
            switch tag
                case 'CyPet'
                    frameInfo(frameNum,9)=obj.ImagePara.CFPexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.CFPexcitationIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,2),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'sfGFP'
                    frameInfo(frameNum,9)=obj.ImagePara.GFPexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.GFPexcitationIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,2),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'sfGFP1'
                    frameInfo(frameNum,9)=obj.ImagePara.GFPexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.GFPexcitationIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,2),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'sfGFP-CyOFP'
                    frameInfo(frameNum,9)=obj.ImagePara.cyOFPexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.cyOFPexcitationIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,2),strcat(feilTag,'\image','sfGFP',num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'Venus'
                    frameInfo(frameNum,9)=obj.ImagePara.YFPexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.YFPexcitationIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,2),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'mAmetrine'
                    frameInfo(frameNum,9)=obj.ImagePara.mAmetrine_exposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.mAmetrine_excitationIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,2),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'mScarletI'
                    frameInfo(frameNum,9)=obj.ImagePara.RFPexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.RFPexcitationIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,1),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'mCherry'
                    frameInfo(frameNum,9)=obj.ImagePara.RFPexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.RFPexcitationIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,1),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'TDmsURFP'
                    frameInfo(frameNum,9)=obj.ImagePara.iRFPexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.iRFPexcitationIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,1),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'CyOFP'
                    frameInfo(frameNum,9)=obj.ImagePara.cyOFPexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.cyOFPexcitationIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,1),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'BlueControl'
                    frameInfo(frameNum,9)=obj.ImagePara.BlueControlexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.BlueControlIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,2),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'GreenControl'
                    frameInfo(frameNum,9)=obj.ImagePara.GreenControlexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.GreenControlIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,2),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'RedControl'
                    frameInfo(frameNum,9)=obj.ImagePara.RedControlexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.RedControlIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,1),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'BF'
                    frameInfo(frameNum,9)=obj.ImagePara.BFexposeTime;
                    frameInfo(frameNum,10)=0;
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageBF(:,:,1),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
                case 'Tracking'
                    if strcmp(obj.xyzC.focalPlan,'BF')
                        imageTracking=obj.currentImage.maskBF(:,:,1);
                    else
                        imageTracking=obj.currentImage.maskFluo(:,:,1);
                    end
                    save(strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.mat'),'imageTracking');
                case 'BlueMask'
                    
                    frameInfo(frameNum,9)=obj.ImagePara.BlueControlexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.BlueControlIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    switch obj.Calib.lightControl
                        case 'dynamic'
                            BlueMask=obj.currentImage.mosaicBlueMask;
                            save(strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.mat'),'BlueMask');
                        case 'static'
                            if   size(frameInfo,1)==1
                                BlueMask=obj.currentImage.mosaicBlueMask;
                                save(strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.mat'),'BlueMask');
                            end
                    end
                    %                     end
                case 'GreenMask'
                    frameInfo(frameNum,9)=obj.ImagePara.GreenControlexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.GreenControlIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    GreenMask=obj.currentImage.mosaicGreenMask;
                    save(strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.mat'),'GreenMask');
                    switch obj.Calib.lightControl
                        case 'dynamic'
                            GreenMask=obj.currentImage.mosaicGreenMask;
                            save(strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.mat'),'GreenMask');
                        case 'static'
                            if   size(frameInfo,1)==1
                                GreenMask=obj.currentImage.mosaicGreenMask;
                                save(strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.mat'),'GreenMask');
                            end
                    end
                case 'RedMask'
                    frameInfo(frameNum,9)=obj.ImagePara.RedControlexposeTime;
                    frameInfo(frameNum,10)=str2double(obj.ImagePara.RedControlIntensity);
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    switch obj.Calib.lightControl
                        case 'dynamic'
                            RedMask=obj.currentImage.mosaicRedMask;
                            save(strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.mat'),'RedMask');
                        case 'static'
                            if   size(frameInfo,1)==1
                                RedMask=obj.currentImage.mosaicRedMask;
                                save(strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.mat'),'RedMask');
                            end
                            %                     end
                        otherwise
                            disp('tag do not exist');
                    end
            end
            
        end
        
        
    end
    
end