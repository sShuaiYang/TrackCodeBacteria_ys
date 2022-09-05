classdef ImageProtocal
    
    properties
        ImagePara=ImageParameter;
        DeviceS=DeviceState;
        Calib=Calibration;
        xyzC=xyzControl;
        currentImage=currentImage;
        imageProcessing=imageProcessing;
        hardware=hardwareSetup;
        feildTag=feildTag;
        imageShowing=imageShowing;
        mmc
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
                obj.mmc = loadMMconfigFile('C:\Program Files\Micro-Manager-1.4.23\allDevices2.cfg');
                disp('Import Fiji..');
                addpath('D:\Fiji\Fiji.app\scripts');
                Miji;
                import ij.IJ;
                disp('Setup hareware...');
                obj.mmc.setConfig('channel','teal');
                obj = obj.getCurrentDeviceState();
                disp('Welcome to a new experiment');
                obj = hardwareSet(obj);
                obj.currentImage=obj.currentImage.initializeCurrentImage();
                obj.imageShowing=obj.imageShowing.initializeFigureSetting();
                obj.Calib=obj.Calib.initializeCalibration();
                obj.experimentIdx=input('Please input experimental idx:');
                obj.experimentType=input('Please input experimental type:1 flow cell, 2 static cell, 3 plate');
                olympus=input('Which light path of olympus will you select:''1x'' or ''1.6x'';');
                lightControl=input('Please input mode of light control:''no'',''static'',''dynamic'';');
                obj.Calib.Olympus=olympus;
                obj.Calib.lightControl=lightControl;
                objective=input('Which objective of olympus will you select:[100];');
                obj.Calib.Objective=objective;
                switch olympus
                    case '1x'
                        obj.Calib.scale=0.065*100/objective;
                    case '1.6x'
                        obj.Calib.scale=0.065/1.6*100/objective;
                end
                if obj.Calib.Objective==20
                    obj.mmc.setProperty('CRISP','Calibration Gain','80');
                    obj.mmc.setProperty('CRISP','GainMultiplier','10');
                    obj.mmc.setProperty('CRISP','Objective NA','0.86');
                end
            catch ME
                obj.mmc.unloadAllDevices();
                error('Initialization failed');
            end
            
            %             obj.mmc=load();                            %
            %             obj.mij=xx;                                %
        end
        
        function obj = hardwareSet(obj)
            obj.mmc.setProperty('Spectra','State','0');
            obj.mmc.setProperty('Wheel-B','Speed','3');
            obj.mmc.setProperty('Wheel-A','Speed','3');
            obj.mmc.setProperty('ZStage','Speed-S','0.3');
            obj.mmc.setProperty('Mosaic3','PixelMode','16GraysLinear');%% '16Grays','16GraysLinear','256Grays','64Grays','BlackAndWhite'
            obj.mmc.setProperty('Mosaic3','OverlapMode','On');%% 'On','Off'
            obj.mmc.setProperty('Mosaic3','GrayScaleExposureTime','0.05');%min exposure time 0.05 s,max exposure time 1 s
            obj.mmc.setProperty('Mosaic3','ExposureTime','100');%%min exposure time 0.05 s
            
            obj.mmc.waitForSystem();
            obj.hardware.ZStageSpeed=char(obj.mmc.getProperty('ZStage','Speed-S'));
            obj.hardware.XYStageSpeed=char(obj.mmc.getProperty('XYStage','Speed-S'));
            obj.hardware.WheelASpeed=char(obj.mmc.getProperty('Wheel-A','Speed'));
            obj.hardware.WheelBSpeed=char(obj.mmc.getProperty('Wheel-B','Speed'));
            obj.hardware.camera1AcquisitionWindow=char(obj.mmc.getProperty('Physical Camera 1','AcquisitionWindow'));
            obj.hardware.camera1DynamicRange=char(obj.mmc.getProperty('Physical Camera 1','Sensitivity/DynamicRange'));
            obj.hardware.camera2AcquisitionWindow=char(obj.mmc.getProperty('Physical Camera 2','AcquisitionWindow'));
            obj.hardware.camera2DynamicRange=char(obj.mmc.getProperty('Physical Camera 2','Sensitivity/DynamicRange'));
            obj.hardware.mosaic3PixelMode=char(obj.mmc.getProperty('Mosaic3','PixelMode'));
            obj.hardware.mosaic3ExposeTime=char(obj.mmc.getProperty('Mosaic3','ExposureTime'));
            obj.hardware.mosaic3GrayScaleExposeTime=char(obj.mmc.getProperty('Mosaic3','GrayScaleExposureTime'));
            obj.hardware.mosaic3OverlapMode=char(obj.mmc.getProperty('Mosaic3','OverlapMode'));
            
            
        end
        
        function obj = getCurrentDeviceState(obj)
            obj.DeviceS.FilterCube = char( obj.mmc.getProperty('FilterCube','Label'));
            state=0;
            while state==0
                try
                    obj.DeviceS.CRISPState = char( obj.mmc.getProperty('CRISP','CRISP State'));
                    state=1;
                catch
                    state=0;
                end
            end
            obj.DeviceS.PhysicalCamera1Trigger = char(obj.mmc.getProperty('Physical Camera 1','TriggerMode'));
            obj.DeviceS.PhysicalCamera1Exposure = char(obj.mmc.getProperty('Physical Camera 1','Exposure'));
            obj.DeviceS.PhysicalCamera2Trigger = char(obj.mmc.getProperty('Physical Camera 2','TriggerMode'));
            obj.DeviceS.PhysicalCamera2Exposure = char(obj.mmc.getProperty('Physical Camera 2','Exposure'));
            obj.DeviceS.CameraIdx = char(obj.mmc.getCameraDevice());
            obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
            obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
            obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
            obj.DeviceS.Mosaic = char(obj.mmc.getSLMDevice());
            obj.DeviceS.LumencorVioletIntensity = char(obj.mmc.getProperty('Spectra','Violet_Level'));
            obj.DeviceS.LumencorBlueIntensity = char(obj.mmc.getProperty('Spectra','Blue_Level'));
            obj.DeviceS.LumencorCyanIntensity = char(obj.mmc.getProperty('Spectra','Cyan_Level'));
            obj.DeviceS.LumencorTealIntensity = char(obj.mmc.getProperty('Spectra','Teal_Level'));
            obj.DeviceS.LumencorGreenIntensity = char(obj.mmc.getProperty('Spectra','Green_Level'));
            obj.DeviceS.LumencorRedIntensity = char(obj.mmc.getProperty('Spectra','Red_Level'));
            %             obj.DeviceS.LumencorVioletEnable = char(obj.mmc.getProperty('Spectra','Violet_Level'));
            %             obj.DeviceS.LumencorBlueEnable = char(obj.mmc.getProperty('Spectra','Blue_Enable'));
            %             obj.DeviceS.LumencorCyanEnable = char(obj.mmc.getProperty('Spectra','Cyan_Enable'));
            %             obj.DeviceS.LumencorTealEnable = char(obj.mmc.getProperty('Spectra','Teal_Enable'));
            %             obj.DeviceS.LumencorGreenEnable = char(obj.mmc.getProperty('Spectra','Green_Enable'));
            %             obj.DeviceS.LumencorRedEnable = char(obj.mmc.getProperty('Spectra','Red_Enable'));
            
            
        end
        
        function obj = setCameraIdx(obj,CameraIdx) %set only first time
            if ~strcmp(obj.DeviceS.CameraIdx,CameraIdx)
                obj.mmc.setCameraDevice(CameraIdx);
                obj.DeviceS.CameraIdx = char(obj.mmc.getCameraDevice());
            end
            switch CameraIdx
                case 'Physical Camera 1'
                    if ~strcmp(char(obj.DeviceS.PhysicalCamera1Trigger),'Internal (Recommended for fast acquisitions)')
                        obj.mmc.setProperty('Physical Camera 1','TriggerMode','Internal (Recommended for fast acquisitions)');
                        obj.DeviceS.PhysicalCamera1Trigger='Internal (Recommended for fast acquisitions)';
                    end
                    if ~strcmp(char(obj.DeviceS.PhysicalCamera2Trigger),'Internal (Recommended for fast acquisitions)')
                        obj.mmc.setProperty('Physical Camera 2','TriggerMode','Internal (Recommended for fast acquisitions)');
                        obj.DeviceS.PhysicalCamera2Trigger='Internal (Recommended for fast acquisitions)';
                    end
                case 'Physical Camera 2'
                    if ~strcmp(char(obj.DeviceS.PhysicalCamera1Trigger),'Internal (Recommended for fast acquisitions)')
                        obj.mmc.setProperty('Physical Camera 1','TriggerMode','Internal (Recommended for fast acquisitions)');
                        obj.DeviceS.PhysicalCamera1Trigger='Internal (Recommended for fast acquisitions)';
                    end
                    if ~strcmp(char(obj.DeviceS.PhysicalCamera2Trigger),'Internal (Recommended for fast acquisitions)')
                        obj.mmc.setProperty('Physical Camera 2','TriggerMode','Internal (Recommended for fast acquisitions)');
                        obj.DeviceS.PhysicalCamera2Trigger='Internal (Recommended for fast acquisitions)';
                    end
                case 'Multi Camera'
                    if ~strcmp(char(obj.DeviceS.PhysicalCamera1Trigger),'Internal (Recommended for fast acquisitions)')
                        obj.mmc.setProperty('Physical Camera 1','TriggerMode','Internal (Recommended for fast acquisitions)');
                        obj.DeviceS.PhysicalCamera1Trigger='Internal (Recommended for fast acquisitions)';
                    end
                    if ~strcmp(char(obj.DeviceS.PhysicalCamera2Trigger),'External')
                        obj.mmc.setProperty('Physical Camera 2','TriggerMode','External');
                        obj.DeviceS.PhysicalCamera2Trigger='External';
                    end
            end
        end
        
        
        function obj=setChannelDeviceState(obj,channel)
            switch channel
                case 'BF'
                    obj.mmc.setProperty('Arduino-Switch','Label','1.ch1');%select BF LED;
                    obj.mmc.setProperty('Arduino-Shutter','OnOff','0');
                    if ~strcmp(obj.DeviceS.FilterCube,'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube = char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(obj.DeviceS.WheelA,'Filter-1')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA =char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(obj.DeviceS.WheelB,'Filter-1')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Arduino-Shutter')
                        obj.mmc.setShutterDevice('Arduino-Shutter');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter =char(obj.mmc.getShutterDevice());
                    end
                    if obj.mmc.getExposure ~=obj.ImagePara.BFexposeTime
                        obj.mmc.setExposure(obj.ImagePara.BFexposeTime);
                    end
                    case 'BFON'
%                     obj.mmc.setProperty('Arduino-Switch','Label','1.ch1');%select BF LED;
%                     obj.mmc.setProperty('Arduino-Shutter','OnOff','0');
                    if ~strcmp(obj.DeviceS.FilterCube,'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube = char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(obj.DeviceS.WheelA,'Filter-1')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA =char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(obj.DeviceS.WheelB,'Filter-1')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'')
                        obj.mmc.setShutterDevice('');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter =char(obj.mmc.getShutterDevice());
                    end
                    if obj.mmc.getExposure ~=obj.ImagePara.BFexposeTime
                        obj.mmc.setExposure(obj.ImagePara.BFexposeTime);
                    end
                case 'VioletFeild'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-1')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-1')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorVioletIntensity,obj.ImagePara.mAmetrine_excitationIntensity)
                        obj.mmc.setProperty('Spectra','Violet_Level',obj.ImagePara.mAmetrine_excitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorBlueIntensity=char(obj.mmc.getProperty('Spectra','Violet_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','1');
                    if obj.mmc.getExposure ~=obj.ImagePara.mAmetrine_exposeTime
                        obj.mmc.setExposure(obj.ImagePara.mAmetrine_exposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'BlueFeild'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-1')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-1')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorBlueIntensity,obj.ImagePara.CFPexcitationIntensity)
                        obj.mmc.setProperty('Spectra','Blue_Level',obj.ImagePara.CFPexcitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorBlueIntensity=char(obj.mmc.getProperty('Spectra','Blue_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','1');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.CFPexposeTime
                        obj.mmc.setExposure(obj.ImagePara.CFPexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'CyanFeild'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-1')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-1')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorCyanIntensity,obj.ImagePara.GFPexcitationIntensity)
                        obj.mmc.setProperty('Spectra','Cyan_Level',obj.ImagePara.GFPexcitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorBlueIntensity=char(obj.mmc.getProperty('Spectra','Cyan_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','1');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.GFPexposeTime
                        obj.mmc.setExposure(obj.ImagePara.GFPexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'TealFeild'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-1')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-1')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorTealIntensity,obj.ImagePara.YFPexcitationIntensity)
                        obj.mmc.setProperty('Spectra','Teal_Level',obj.ImagePara.YFPexcitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorBlueIntensity=char(obj.mmc.getProperty('Spectra','Teal_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','1');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.YFPexposeTime
                        obj.mmc.setExposure(obj.ImagePara.YFPexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'GreenFeild'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-1')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-1')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorGreenIntensity,obj.ImagePara.RFPexcitationIntensity)
                        obj.mmc.setProperty('Spectra','Green_Level',obj.ImagePara.RFPexcitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorBlueIntensity=char(obj.mmc.getProperty('Spectra','Green_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','1');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.RFPexposeTime
                        obj.mmc.setExposure(obj.ImagePara.RFPexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'RedFeild'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-1')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-1')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorRedIntensity,obj.ImagePara.iRFPexcitationIntensity)
                        obj.mmc.setProperty('Spectra','Red_Level',obj.ImagePara.iRFPexcitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorBlueIntensity=char(obj.mmc.getProperty('Spectra','Red_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','1');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.iRFPexposeTime
                        obj.mmc.setExposure(obj.ImagePara.iRFPexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'GFPCyOFP'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-2')
                        obj.mmc.setProperty('FilterCube','Label','Position-2');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-9')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-9');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-0')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-0');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorCyanIntensity,obj.ImagePara.cyOFPexcitationIntensity)
                        obj.mmc.setProperty('Spectra','Cyan_Level',obj.ImagePara.cyOFPexcitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorCyanIntensity=char(obj.mmc.getProperty('Spectra','Cyan_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','1');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.cyOFPexposeTime
                        obj.mmc.setExposure(obj.ImagePara.cyOFPexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'PVD'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-2')
                        obj.mmc.setProperty('FilterCube','Label','Position-2');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-2')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-2');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorVioletIntensity,obj.ImagePara.mAmetrine_excitationIntensity)
                        obj.mmc.setProperty('Spectra','Violet_Level',obj.ImagePara.mAmetrine_excitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorRedIntensity=char(obj.mmc.getProperty('Spectra','Violet_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','1');
                    if obj.mmc.getExposure ~=obj.ImagePara.mAmetrine_exposeTime
                        obj.mmc.setExposure(obj.ImagePara.mAmetrine_exposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'mAmetrine'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-2')
                        obj.mmc.setProperty('FilterCube','Label','Position-2');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-9')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-9');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorVioletIntensity,obj.ImagePara.mAmetrine_excitationIntensity)
                        obj.mmc.setProperty('Spectra','Violet_Level',obj.ImagePara.mAmetrine_excitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorRedIntensity=char(obj.mmc.getProperty('Spectra','Violet_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','1');
                    if obj.mmc.getExposure ~=obj.ImagePara.mAmetrine_exposeTime
                        obj.mmc.setExposure(obj.ImagePara.mAmetrine_exposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'GFP'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-2')
                        obj.mmc.setProperty('FilterCube','Label','Position-2');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-9')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-9');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-9')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-9');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorCyanIntensity,obj.ImagePara.GFPexcitationIntensity)
                        obj.mmc.setProperty('Spectra','Cyan_Level',obj.ImagePara.GFPexcitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorCyanIntensity=char(obj.mmc.getProperty('Spectra','Cyan_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','1');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.GFPexposeTime
                        obj.mmc.setExposure(obj.ImagePara.GFPexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'CyOFP'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-2')
                        obj.mmc.setProperty('FilterCube','Label','Position-2');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-0')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-0');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorCyanIntensity,obj.ImagePara.cyOFPexcitationIntensity)
                        obj.mmc.setProperty('Spectra','Cyan_Level',obj.ImagePara.cyOFPexcitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorCyanIntensity=char(obj.mmc.getProperty('Spectra','Cyan_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','1');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.cyOFPexposeTime
                        obj.mmc.setExposure(obj.ImagePara.cyOFPexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'RFP'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-2')
                        obj.mmc.setProperty('FilterCube','Label','Position-2');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-9')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-9');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorGreenIntensity,obj.ImagePara.RFPexcitationIntensity)
                        obj.mmc.setProperty('Spectra','Green_Level',obj.ImagePara.RFPexcitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorGreenIntensity=char(obj.mmc.getProperty('Spectra','Green_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','1');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.RFPexposeTime
                        obj.mmc.setExposure(obj.ImagePara.RFPexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'iRFP'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-2')
                        obj.mmc.setProperty('FilterCube','Label','Position-2');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-8')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-8');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorRedIntensity,obj.ImagePara.iRFPexcitationIntensity)
                        obj.mmc.setProperty('Spectra','Red_Level',obj.ImagePara.iRFPexcitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorRedIntensity=char(obj.mmc.getProperty('Spectra','Red_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','1');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.iRFPexposeTime
                        obj.mmc.setExposure(obj.ImagePara.iRFPexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'CFP'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-1')
                        obj.mmc.setProperty('FilterCube','Label','Position-1');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-0')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-0');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorBlueIntensity,obj.ImagePara.CFPexcitationIntensity)
                        obj.mmc.setProperty('Spectra','Blue_Level',obj.ImagePara.CFPexcitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorRedIntensity=char(obj.mmc.getProperty('Spectra','Blue_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','1');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.CFPexposeTime
                        obj.mmc.setExposure(obj.ImagePara.CFPexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'YFP'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-1')
                        obj.mmc.setProperty('FilterCube','Label','Position-1');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-7')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-7');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorTealIntensity,obj.ImagePara.YFPexcitationIntensity)
                        obj.mmc.setProperty('Spectra','Teal_Level',obj.ImagePara.YFPexcitationIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorRedIntensity=char(obj.mmc.getProperty('Spectra','Teal_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','1');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.YFPexposeTime
                        obj.mmc.setExposure(obj.ImagePara.YFPexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicWhiteImage)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicWhiteImage);
                        MIJ.run('Close All');
                    end
                case 'BlueControl'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-1')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-1')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorBlueIntensity,obj.ImagePara.BlueControlIntensity)
                        obj.mmc.setProperty('Spectra','Blue_Level',obj.ImagePara.BlueControlIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorBlueIntensity=char(obj.mmc.getProperty('Spectra','Blue_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','1');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.BlueControlexposeTime
                        obj.mmc.setExposure(obj.ImagePara.BlueControlexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicBlueMask)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicBlueMask);
                        MIJ.run('Close All');
                    end
                case 'BlueControlTest'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-9')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-9');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-0')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-0');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.mmc.getProperty('Spectra','Blue_Level'),obj.ImagePara.BlueControlIntensity)
                        obj.mmc.setProperty('Spectra','Blue_Level',obj.ImagePara.BlueControlIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorBlueIntensity=char(obj.mmc.getProperty('Spectra','Blue_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','1');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.BlueControlexposeTime
                        obj.mmc.setExposure(obj.ImagePara.BlueControlexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicBlueMask)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicBlueMask);
                        MIJ.run('Close All');
                    end
                    
                case 'GreenControl'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-1')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-1')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorTealIntensity,obj.ImagePara.GreenControlIntensity)
                        obj.mmc.setProperty('Spectra','Teal_Level',obj.ImagePara.GreenControlIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorGreenIntensity=char(obj.mmc.getProperty('Spectra','Teal_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','1');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.GreenControlexposeTime
                        obj.mmc.setExposure(obj.ImagePara.GreenControlexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicGreenMask)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicGreenMask);
                    end
                case 'GreenControlTest'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-9')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-9');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-9')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-9');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorTealIntensity,obj.ImagePara.GreenControlIntensity)
                        obj.mmc.setProperty('Spectra','Teal_Level',obj.ImagePara.GreenControlIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorGreenIntensity=char(obj.mmc.getProperty('Spectra','Teal_Level'));
                    end
                    obj.mmc.setProperty('Spectra','Teal_Enable','1');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','0');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    if obj.mmc.getExposure ~=obj.ImagePara.GreenControlexposeTime
                        obj.mmc.setExposure(obj.ImagePara.GreenControlexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicGreenMask)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicGreenMask);
                    end
                case 'RedControl'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-1')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-1')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorRedIntensity,obj.ImagePara.RedControlIntensity)
                        obj.mmc.setProperty('Spectra','Red_Level',obj.ImagePara.RedControlIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorRedIntensity=char(obj.mmc.getProperty('Spectra','Red_Level'));
                    end
                    
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','1');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    
                    if obj.mmc.getExposure ~=obj.ImagePara.RedControlexposeTime
                        obj.mmc.setExposure(obj.ImagePara.RedControlexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicRedMask)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicRedMask);
                    end
                case 'RedControlTest'
                    if ~strcmp(char(obj.DeviceS.FilterCube),'Position-3')
                        obj.mmc.setProperty('FilterCube','Label','Position-3');
                        obj.mmc.waitForDevice('FilterCube');
                        obj.DeviceS.FilterCube =char(obj.mmc.getProperty('FilterCube','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelA),'Filter-1')
                        obj.mmc.setProperty('Wheel-A','Label','Filter-1');
                        obj.mmc.waitForDevice('Wheel-A');
                        obj.DeviceS.WheelA = char(obj.mmc.getProperty('Wheel-A','Label'));
                    end
                    if ~strcmp(char(obj.DeviceS.WheelB),'Filter-0')
                        obj.mmc.setProperty('Wheel-B','Label','Filter-0');
                        obj.mmc.waitForDevice('Wheel-B');
                        obj.DeviceS.WheelB = char(obj.mmc.getProperty('Wheel-B','Label'));
                    end
                    if ~strcmp(obj.DeviceS.Shutter,'Spectra')
                        obj.mmc.setShutterDevice('Spectra');
                        obj.mmc.waitForDevice('Core');
                        obj.DeviceS.Shutter = char(obj.mmc.getShutterDevice());
                    end
                    if ~strcmp(obj.DeviceS.LumencorRedIntensity,obj.ImagePara.RedControlIntensity)
                        obj.mmc.setProperty('Spectra','Red_Level',obj.ImagePara.RedControlIntensity);
                        obj.mmc.waitForDevice('Spectra');
                        obj.DeviceS.LumencorRedIntensity=char(obj.mmc.getProperty('Spectra','Red_Level'));
                    end
                    
                    obj.mmc.setProperty('Spectra','Teal_Enable','0');
                    obj.mmc.setProperty('Spectra','Green_Enable','0');
                    obj.mmc.setProperty('Spectra','Blue_Enable','0');
                    obj.mmc.setProperty('Spectra','Cyan_Enable','0');
                    obj.mmc.setProperty('Spectra','Red_Enable','1');
                    obj.mmc.setProperty('Spectra','Violet_Enable','0');
                    
                    if obj.mmc.getExposure ~=obj.ImagePara.RedControlexposeTime
                        obj.mmc.setExposure(obj.ImagePara.RedControlexposeTime);
                    end
                    if ~isempty(obj.currentImage.mosaicRedMask)
                        setImage2Mosaic(obj.mmc,obj.currentImage.mosaicRedMask);
                    end
                    
                    
            end
        end
        
        function obj= getSnapImage(obj,channelIdx)
            switch obj.DeviceS.CameraIdx
                case 'Physical Camera 1'
                    
                case 'Physical Camera 2'
                    
                    
                case 'Multi Camera'
                    switch channelIdx
                        case 'BF'
                            %                     obj.mmc.setProperty('Arduino-Shutter','OnOff','1');%BF LED on
                            obj.mmc.snapImage(); %after set shutter to Arduino-Switch,run mmc.snapImage will tigger BF LED on
                            obj.mmc.setProperty('Arduino-Shutter','OnOff','0');%BF LED off, setting Arduino-Shutter off will faster turnoff the BF LED
                        case 'Mosaic'
                            obj.mmc.displaySLMImage(obj.DeviceS.Mosaic);
                            obj.mmc.snapImage();
                        case 'BFON'
                            %                     obj.mmc.setProperty('Arduino-Shutter','OnOff','1');%BF LED on
                            obj.mmc.snapImage(); %after set shutter to Arduino-Switch,run mmc.snapImage will tigger BF LED on
%                             obj.mmc.setProperty('Arduino-Shutter','OnOff','0');%BF LED off, setting Arduino-Shutter off will faster turnoff the BF LED    
                        otherwise
                            disp('None channelIdx');
                    end
                    
                    image1 = obj.mmc.getImage(0);  % Physical Camera 1 get image
                    image2 = obj.mmc.getImage(1);  % Physical Camera 2 get image
                    width = obj.mmc.getImageWidth();
                    height = obj.mmc.getImageHeight();
                    if obj.mmc.getBytesPerPixel == 2
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
                        if ~isempty(obj.currentImage.imageBFOri);
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
                        if ~isempty(obj.currentImage.imageFluoOri);
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
        
        function obj = moveFocalPlane(obj,focalPlane,isAutoFocus)
            zMovePrecision=0.1; %0.1 um
            if strcmp(obj.xyzC.basePlan,focalPlane)
                state=0;
                while state==0
                    try
                        obj.mmc.setProperty('CRISP','CRISP State','Lock');
                        state=1;
                    catch
                        state=0;
                    end
                end
                return;
            end
            state=0;
            while state==0
                try
                    obj.DeviceS.CRISPState= char( obj.mmc.getProperty('CRISP','CRISP State'));
                    state=1;
                catch
                    state=0;
                end
            end
            state=0;
            if ~strcmp(focalPlane,obj.xyzC.focalPlan)
                if strcmp(obj.DeviceS.CRISPState,'In Focus')||strcmp(obj.DeviceS.CRISPState,'Lock')
                    while state==0
                        try
                            obj.mmc.setProperty('CRISP','CRISP State','Ready');
                            obj.DeviceS.CRISPState = obj.mmc.getProperty('CRISP','CRISP State');
                            state=1;
                        catch
                            state=0;
                        end
                    end
                end
                
                if strcmp(obj.xyzC.focalPlan,'BF')
                    zshift=obj.xyzC.BFzShift;
                elseif strcmp(obj.xyzC.focalPlan,'FL')
                    zshift=obj.xyzC.FLzShift;
                elseif strcmp(obj.xyzC.focalPlan,'BlueControl')
                    zshift=obj.Calib.BlueZshift;
                elseif strcmp(obj.xyzC.focalPlan,'GreenControl')
                    zshift=obj.Calib.GreenZshift;
                elseif strcmp(obj.xyzC.focalPlan,'RedControl')
                    zshift=obj.Calib.RedZshift;
                end
                
                switch focalPlane
                    case 'BF'
                        if abs(obj.xyzC.BFzShift-zshift)>zMovePrecision
                            zposition=obj.mmc.getPosition();
                            obj.setZPostion(obj.xyzC.BFzShift-zshift);
                            obj.xyzC.focalPlan='BF';
                        end
                    case 'FL'
                        if abs(obj.xyzC.FLzShift-zshift)>zMovePrecision
                            zposition=obj.mmc.getPosition();
                            obj.setZPostion(obj.xyzC.FLzShift-zshift);
                            obj.xyzC.focalPlan='FL';
                        end
                    case 'BlueControl'
                        if abs(obj.Calib.BlueZshift-zshift)>zMovePrecision
                            zposition=obj.mmc.getPosition();
                            obj.setZPostion(obj.Calib.BlueZshift-zshift);
                            obj.xyzC.focalPlan='BlueControl';
                        end
                    case 'GreenControl'
                        if abs(obj.Calib.GreenZshift-zshift)>zMovePrecision
                            zposition=obj.mmc.getPosition();
                            obj.setZPostion(obj.Calib.GreenZshift-zshift);
                            obj.xyzC.focalPlan = 'GreenControl';
                        end
                    case 'RedControl'
                        if abs(obj.Calib.RedZshift-zshift)>zMovePrecision
                            zposition=obj.mmc.getPosition();
                            obj.setZPostion(obj.Calib.RedZshift-zshift);
                            obj.xyzC.focalPlan = 'RedControl';
                        end
                end
                
                if isAutoFocus
                    obj.mmc.setProperty('CRISP','CRISP State','Idle');
                    pause(0.5);
                    obj.mmc.setProperty('CRISP','CRISP State','loG_cal');
                    pause(3);
                    obj.mmc.setProperty('CRISP','CRISP State','gain_Cal');
                    pause(3);
                    obj.mmc.setProperty('CRISP','CRISP State','Lock');
                    pause(2);
                    state=0;
                    while state==0
                        try
                            obj.DeviceS.CRISPState = obj.mmc.getProperty('CRISP','CRISP State');
                            state=1;
                        catch
                            state=0;
                        end
                    end
                    obj.xyzC.focalPlan=focalPlane;
                else
                    obj.xyzC.focalPlan='BF';
                end
            end
        end
        
        
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
            calibFeildNum=size(obj.Calib.calibrationFeild,1);
            obj.Calib.calibrationFeild(calibFeildNum+1,1)=obj.mmc.getXPosition();
            obj.Calib.calibrationFeild(calibFeildNum+1,2)=obj.mmc.getYPosition();
            obj=obj.getCurrentDeviceState();
            if ~isempty(lightSelect)
                for i=1:numel(lightSelect);
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
            
            if isempty(obj.feildTag.name)
                obj.feildTag.name{1}='feild0001';
                obj.feildTag.xyPosition(1,1)=obj.mmc.getXPosition();
                obj.feildTag.xyPosition(1,2)=obj.mmc.getYPosition();
                obj.feildTag.tag{1}='none';
                obj.feildTag.tagValue(1,:)=0;
            else
                currentFeildTag=size(obj.feildTag.name,2);
                obj.feildTag.name{currentFeildTag+1}=strcat('feild',num2str(currentFeildTag+1,'%.4d'));
                obj.feildTag.xyPosition(currentFeildTag+1,1)=obj.mmc.getXPosition();
                obj.feildTag.xyPosition(currentFeildTag+1,2)=obj.mmc.getYPosition();
                obj.feildTag.tag{currentFeildTag+1}='none';
                obj.feildTag.tagValue(currentFeildTag+1,:)=0;
            end
        end
        
        function obj=addMontageFeild(obj,xnumber,ynumber,gap)
            %             if ~isempty(obj.feildTag.name)
            %                 obj.feildTag.name={};
            %                 obj.feildTag.xyPosition=[];
            %                 obj.feildTag.tag={};
            %                 obj.feildTag.tagValue=[];
            %             end
            currentFeildTag1=size(obj.feildTag.name,2);
            imageSize=[2048,2048];
            currentXYPosition(1)=obj.mmc.getXPosition();
            currentXYPosition(2)=obj.mmc.getYPosition();
            for x=1:xnumber
                for y=1:ynumber
                    
                    currentFeildTag=(x-1)*ynumber+y+currentFeildTag1;
                    obj.feildTag.name{currentFeildTag}=strcat('feild',num2str(currentFeildTag,'%04d'));
                    obj.feildTag.xyPosition(currentFeildTag,1:2)=[currentXYPosition(1)+(1+gap)*(x-1)*imageSize(1)*obj.Calib.scale,currentXYPosition(2)+(1+gap)*(y-1)*imageSize(2)*obj.Calib.scale];
                    obj.feildTag.tag{currentFeildTag}='none';
                    obj.feildTag.tagValue(currentFeildTag,:)=0;
                end
            end
            
            
        end
        
        function obj=deleteXYPosition(obj,varargin)
            if isempty(varargin)
                obj.feildTag.name={};
                obj.feildTag.xyPosition=[];
                obj.feildTag.tag={};
                obj.feildTag.tagValue=[];
            else
                delectXYFeild=varargin{1};
                delectXYFeild=sort(delectXYFeild,'descend');
                for iFeild=1:numel(delectXYFeild)
                    obj.feildTag.name(delectXYFeild(iFeild))=[];
                    obj.feildTag.xyPosition(delectXYFeild(iFeild),:)=[];
                    obj.feildTag.tag(delectXYFeild(iFeild))=[];
                    obj.feildTag.tagValue(delectXYFeild(iFeild))=[];
                    if delectXYFeild(iFeild)<size(obj.feildTag.name,1)+1
                        for i=delectXYFeild(iFeild):size(obj.feildTag.name,2)
                            obj.feildTag.name{i}=strcat('feild',num2str(i,'%.4d'));
                        end
                    end
                    
                end
            end
        end
        
        function obj=showFeildPosition(obj)
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
            if ~isempty(obj.feildTag.name)
                for i=1:size(obj.feildTag.xyPosition,1)
                    label=num2str(i);
                    plot(obj.feildTag.xyPosition(i,1),obj.feildTag.xyPosition(i,2));
                    text(obj.feildTag.xyPosition(i,1),obj.feildTag.xyPosition(i,2),label,'Color',[1,0,0]);
                    rectangle('Position',[obj.feildTag.xyPosition(i,1)-rectangleW/2,obj.feildTag.xyPosition(i,2)-rectangleH/2,rectangleW,rectangleH], 'LineWidth',2,...
                        'LineStyle','--',...
                        'EdgeColor',[0 0.498039215803146 0]);
                end
            end
            
        end
        
        function obj=move2Feild(obj,feildNum)
            obj.mmc.setXYPosition(obj.feildTag.xyPosition(feildNum,1),obj.feildTag.xyPosition(feildNum,2));
            %             obj.mmc.waitForDevice('XYStage');
            obj.xyzC.currentFeild=obj.feildTag.name{feildNum};
            obj.xyzC.XYpositions(1)=obj.mmc.getXPosition();
            obj.xyzC.XYpositions(2)=obj.mmc.getYPosition();
        end
        
        function obj=saveImageFile(obj,tag,imageProcessing)
            
            currentFeild=obj.xyzC.currentFeild;
            currentFeildFile=strcat(obj.myDir,'\',currentFeild);
            if ~exist(currentFeildFile,'file')
                mkdir(currentFeildFile);
            end
            feilTag=strcat(currentFeildFile,'\',tag);
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
            frameInfo(frameNum,7)=obj.mmc.getXPosition();
            frameInfo(frameNum,8)=obj.mmc.getYPosition();
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
                     case 'BFOri'
                    frameInfo(frameNum,9)=obj.ImagePara.BFexposeTime;
                    frameInfo(frameNum,10)=0;
                    save(strcat(feilTag,'\frameInfo.mat'),'frameInfo');
                    imwrite(obj.currentImage.imageOri(:,:,2),strcat(feilTag,'\image',tag,num2str(frameNum,'%05d'),'.tif'),'tif');
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
