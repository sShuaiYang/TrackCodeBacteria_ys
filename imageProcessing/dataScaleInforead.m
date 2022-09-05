function dataOut=dataScaleInforead(filein)
% input filein is the URL of the information.txt
fidin=fopen(filein,'r');
nline=0;
while ~feof(fidin) % 判断是否为文件末尾
    tline=fgetl(fidin); % 从文件读行
    nline=nline+1;
    if nline==5
        % find the line *** x : 1721 * 0.061728 : um ***
        info=textscan(tline,'%s%s%f%s%f');
        scaleInfo=info{5};
        xSize=info{3};
    end
    if nline==6
        info=textscan(tline,'%s%s%f%s%f');
        ySize=info{3};
    end
    %find the line *** Repeat T - 20000 times (5 sec) ***
    isRepeatLine=strfind(tline,'Repeat');
    if ~isempty(isRepeatLine)
        info=textscan(tline,'%s %s %s %f %s %s %s');
        if strcmp(info{7}{1}(1:end-1),'sec')
            timeInterval=str2num(info{6}{1}(2:end));
        end
        if strcmp(info{7}{1}(1:end-1),'ms')
           timeInterval=str2num(info{6}{1}(2:end))/1000;
        end
    end
end
fclose(fidin);
dataOut.rcSize=[ySize,xSize];
dataOut.scaleInfo=scaleInfo;
dataOut.timeInterval=timeInterval;
end