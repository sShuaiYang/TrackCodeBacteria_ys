function BleRGBSwitch_CTControl(LEDIntesntiy_96wells, CT, dutyCycle)
%% Wells96 RGB 光不同强度、不同周期CT和占空比dutyCycle控制
% LEDIntesntiy_96wells 8*12 矩阵:范围[0-255];
% CT 8*12 矩阵；dutyCycle  8*12 矩阵；
% 每一个well的周期和占空比都可以自定义
% CT = repmat((1:12)*10,8,1); 每一列同样周期
% 常开或者常关的周期暂时定为数值1
% dutyCycle = repmat(ones(1,12)*0.05,8,1);
% dutyCycle(:,end) = 1;dutyCycle(:,1) = 0;
% 占空比0.05；0 代表长关;1代表常开;
% 参考lightControl_CTvariation %Nilei code
% *Shuai Yang 2022/6/14* 
%% Wells 96 LED灯是反转扣下去的 矩阵需要转换
t_start = clock; 
save([datestr(t_start,'yyyymmddTHHMMSS'),'_','t_start.mat'],'t_start');
save([datestr(t_start,'yyyymmddTHHMMSS'),'_','CT.mat'],'CT');
save([datestr(t_start,'yyyymmddTHHMMSS'),'_','dutyCycle.mat'],'dutyCycle');
save([datestr(t_start,'yyyymmddTHHMMSS'),'_','LEDIntesntiy_96wells.mat'],'LEDIntesntiy_96wells');

LEDIntesntiy_96wells = fliplr(LEDIntesntiy_96wells);
CT = fliplr(CT);
dutyCycle = fliplr(dutyCycle);

%%  蓝牙连接
clear BleHandle
BleHandle = ble("DC46751FA97E"); % 蓝牙装置的Address
TXInterface = characteristic(BleHandle, ...
    "8653000A-43E6-47B7-9CB0-5FC21D4AE340","8653000C-43E6-47B7-9CB0-5FC21D4AE340");
% 计算每个孔的光照时间余数
onTime = dutyCycle.*CT;

%%  根据Time调用函数计算亮度

tic;
for i= 1:1000000 
    tRun = toc  % Elapsed time
    intensityMat = uint8(zeros(8,12));
    TF = onTime>mod(tRun,CT);
    intensityMat(TF) = LEDIntesntiy_96wells(TF);
    intensityMat = intensityMat(:)'; %intensityMat 需要是8bit 行向量
    data = [0x03 intensityMat 0x0d];% 'b'
    % data = [0x01 intensityMat' 13];% 'g'
    % data = [0x02 intensityMat' 13];% 'r'
    write(TXInterface,data);      
    pause(0.5);
    t_end = clock;
%     if etime(t_end,t_start) > 1440 %控制照射时间
%         break
%     end
end

%%  亮度回0
intensityMat = 0*intensityMat;
data = [0x03 intensityMat 0x0d];
write(TXInterface,data);
disp('back to off')

save([datestr(t_end,'yyyymmddTHHMMSS'),'_','t_end.mat'],'t_end')
end