function dataOut=scaleInforeadMontage(filein)
% input filein is the URL of the information.txt
fidin=fopen(filein,'r');
nline=0;
while ~feof(fidin) % 判断是否为文件末尾
    tline=fgetl(fidin); % 从文件读行
    nline=nline+1;
    switch nline
        case 5
            % find the line *** x : 1721 * 0.061728 : um ***
            info=textscan(tline,'%s%s%f%s%f');
            scaleInfo=info{5};
            xSize=info{3};
        case 6
            info=textscan(tline,'%s%s%f%s%f');
            ySize=info{3};
        case 7
            info=textscan(tline,'%s%s%f%s%f');
            WavelengthNum = info{3};
        case 8
            info=textscan(tline,'%s%s%f%s%f');
            MontageNum = info{3};
        case 9
            info=textscan(tline,'%s%s%f%s%f');
            TimeNum = info{3};
    end
    %find the line *** Repeat T - 20000 times (5 sec) ***
    isRepeatLine=strfind(tline,'Repeat');
    if ~isempty(isRepeatLine)
        info=textscan(tline,'%s %s %s %f %s %s %s');
        if ~isempty(info{7})
        switch info{7}{1}(1:end-1)
            case 'sec'
                timeInterval=str2num(info{6}{1}(2:end));                
            case 'ms'
                timeInterval=str2num(info{6}{1}(2:end))/1000;
            case 'min'
                timeInterval=str2num(info{6}{1}(2:end))*60;
            case 'h'
                timeInterval=str2num(info{6}{1}(2:end))*3600;
        end
        end
    end
end
timeInterval=1;
fclose(fidin);
dataOut.rcSize=[ySize,xSize];
dataOut.scaleInfo=scaleInfo;
dataOut.timeInterval=timeInterval;
dataOut.WavelengthNum=WavelengthNum;
dataOut.MontageNum=MontageNum;
dataOut.TimeNum=TimeNum;
end
