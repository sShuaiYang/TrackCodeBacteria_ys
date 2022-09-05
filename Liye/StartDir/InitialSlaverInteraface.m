function HOSTI=InitialSlaverInteraface(HostIP,PortNum)
disp('helloworld!');
HOSTI=tcpip(HostIP,PortNum, 'NetworkRole', 'client');
fopen(HOSTI);
disp('Connected with HostIP!')
end