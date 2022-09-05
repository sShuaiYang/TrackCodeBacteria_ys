HostIP='192.168.1.12';  
TargetIP='192.168.1.12';
TargetPort=30000;

%%
MessageI=tcpip(TargetIP, 9500 , 'NetworkRole', 'client');
fopen(MessageI);
Message=['MATLAB -r cd(''C:\StartDir'');HostI=SlaverInteraface(''' HostIP ''',' num2str(TargetPort) ');HostI=Standby(HostI);StartListen; &'];
fwrite(MessageI,Message);
fclose(MessageI);
SlaverI=tcpip('0.0.0.0',TargetPort, 'NetworkRole', 'server');
fopen(SlaverI);

MesStruct.Method='BFprocess';
MesStruct.FileDir='E:\2019-12-26-deltpslpelfliC_IP31 ys';
MesStruct.TotalFrame=10;
MesStruct.SaveType='.tif';
MesString=jsonencode(MesStruct);
fwrite(SlaverI,MesString);
