classdef SlaverInteraface
    
    properties
        TCPIPInterface;
        FilePath;
        TaskStruct;
        FieldInfo;
        FileToPInfo;
        FilePedInfo;
        TotalImageNumber;
        ProcessedImageNumber;
    end
    
    methods
        function obj = SlaverInteraface(HostIP,PortNum)
            %UNTITLED 构造此类的实例
            %   此处显示详细说明
            obj.TCPIPInterface=tcpip(HostIP,PortNum, 'NetworkRole', 'client');
            try
                fopen(obj.TCPIPInterface);
                disp(['Connected with host  ' HostIP ':' num2str(PortNum)]);
                parfor i=1:32
                    A=1;
                end
                fwrite(obj.TCPIPInterface,'OK');
            catch ME
                disp(['Failed to connect with host ' HostIP ':' num2str(PortNum)]);
            end
        end
        function disp(obj)
            disp(obj.TCPIPInterface);
        end
        
        function obj=Standby(obj)
            disp('Waiting for task.')
            try
                while(1)
                    if obj.TCPIPInterface.BytesAvailable
                        pause(0.02);
                        Message=char(fread(obj.TCPIPInterface,obj.TCPIPInterface.BytesAvailable)');
                        break
                    end
                end
            catch ME
                disp(ME);
            end
            MethodStruct=jsondecode(Message);
            switch MethodStruct.Method
                case 'BFprocess'
                    MethodStruct
                    obj.TaskStruct=MethodStruct;
                    obj=BuildTask(obj);
            end
        end
        function obj=BuildTask(obj)
            obj.FilePath=obj.TaskStruct.FileDir;
            FileList=dir([obj.FilePath,'\field*']);
            for iField=1:size(FileList,1)
                obj.FieldInfo(iField).FieldNum=str2double(FileList(iField).name(6:end));
                obj.FieldInfo(iField).FieldDir=[FileList(iField).folder '\' FileList(iField).name '\'];
                obj.FieldInfo(iField).TotalFrame=obj.TaskStruct.TotalFrame;
                obj.FieldInfo(iField).ProcessedFrameNum=0;
                obj.FieldInfo(iField).ProcessedFrameInfo=[];
            end
            obj.TotalImageNumber=obj.TaskStruct.TotalFrame*size(FileList,1);
            obj.ProcessedImageNumber=0;
        end
        function obj=SearchAndProcess(obj)
            obj.FileToPInfo=[];
            for iField=1:size(obj.FieldInfo,2)
                FileList=dir([obj.FieldInfo(iField).FieldDir 'BF\*' obj.TaskStruct.SaveType]);
                for iFile=1:size(FileList,1)
                    TempFileStruct.FieldIndex=iField;
                    TempFileStruct.FileName=[FileList(iFile).folder '\' FileList(iFile).name '\'];
                    TempFileStruct.FileNum=str2double(FileList(iFile).name(end-8:end-4));
                    if TempFileStruct.FileNum>obj.FieldInfo(iField).ProcessedFrameNum
                        if isempty(obj.FileToPInfo)
                            obj.FileToPInfo=TempFileStruct;
                        else
                            obj.FileToPInfo(end+1)=TempFileStruct;
                        end
                    end
                end
            end
            for iFile=1:size(obj.FileToPInfo,2)
                if obj.FileToPInfo(iFile).FileNum==1
                    %%  对于第一张的处理
                else
                    %%  对于其他张的处理
                end
                if isempty(obj.FilePedInfo)
                    obj.FilePedInfo=obj.FileToPInfo(iFile);
                else
                    obj.FilePedInfo(end+1)=obj.FileToPInfo(iFile);
                end
            end
            obj.ProcessedImageNumber=size(obj.FilePedInfo,2);
        end
        
    end
end


