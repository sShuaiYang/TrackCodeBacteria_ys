function lightControl_CTvariation(CT,dutyCycle)
%% CT赋值采取12列式的赋值即（0.1，0.1，0.1，0.1，0.1，0.1，0.1，0.1，0.1，0.1，0.1，0.1）；
%% 每一列的占空比8个平行的相等
% Nilei code 适用于每一列同样的周期和占空比控制 % Shuai Yang 
%%  蓝牙连接

% dutyCycle = 0.05;%%固定占空比
clear BleHandle
% BleHandle=ble("DD3BA10790F1");
BleHandle=ble("DC46751FA97E");

TXInterface=characteristic(BleHandle,"8653000A-43E6-47B7-9CB0-5FC21D4AE340","8653000C-43E6-47B7-9CB0-5FC21D4AE340");

%%计算每个孔的光照时间余数
onTime = zeros(1,96);
for i=1:length(CT)
    onTime((i-1)*8+1:i*8)=dutyCycle(i)*CT(i);
end


%%   根据Time调用IntensityFunction函数计算亮度
tic;
for i= 1:1000000
Time = toc
intensityMat = dutyCycleFunction(Time,onTime,CT);
data=[0x03 intensityMat 0x0d];
write(TXInterface,data);
pause(0.5);
end

%%  亮度回0
intensityMat=0*intensityMat;
data=[0x03 intensityMat 0x0d];
write(TXInterface,data);
disp('back to off')
end

function intensityMat=dutyCycleFunction(Time,onTime,CT)
intensityMat = uint8(zeros(1,96));
for i=1:length(CT)
    tempOnTime = onTime((i-1)*8+1:i*8);
tempLogic = tempOnTime>mod(Time,CT(i));
tempIntMat = intensityMat((i-1)*8+1:i*8);
tempIntMat(tempLogic)=20;
intensityMat((i-1)*8+1:i*8) = tempIntMat;
end
end