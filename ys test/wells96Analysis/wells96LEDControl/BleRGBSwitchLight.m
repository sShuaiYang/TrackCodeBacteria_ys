function BleRGBSwitchLight(boardID,color,intensity)
% intensity 范围0~255, 长度为1或96，长度为1时为全部设置为统一亮度，为矩阵时应为8x12
%  示例： BleRGBSwitchLight('C008F215A615','r',0)
%         BleRGBSwitchLight('C008F215A615','g',1:96)
blehandle=ble(boardID);
if length(intensity)==1
    intensitystring=uint8(ones(96,1)*intensity);
else
    intensitystring=uint8(intensity(:));
end
writehandle = blehandle.characteristic('8653000A-43E6-47B7-9CB0-5FC21D4AE340','8653000C-43E6-47B7-9CB0-5FC21D4AE340');
if strcmpi(color,'r')
messtring=[0x02 intensitystring' 13];
end
if strcmpi(color,'g')
messtring=[0x01 intensitystring' 13];
end
if strcmpi(color,'b')
messtring=[0x03 intensitystring' 13];% whole 
% messtring=[0x14 intensitystring' 13];%high frequency 1ms:9ms
% messtring=[0x15 intensitystring' 13];%mid frequency 2ms :18ms
% messtring=[0x16 intensitystring' 13];%low frequency 3ms:27ms
end
write(writehandle,messtring);
clear blehandle writehandle
end