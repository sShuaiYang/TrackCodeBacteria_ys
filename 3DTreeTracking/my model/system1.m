function result=system1()
%　定义的是这样一个自我分解的系统A，规则如下
%  每一段时间，会生产出一定量的A （p0 个/timeInterval）
%  同时A自己有一定的概率会变坏，变成A_(丧失功能) （p3 /timeInterval）
%  A可以分解坏的A_  A+A_=A ,
%  假设在这段时间内，两个蛋白相遇的概率为p1(/timeInterval),那么遇到坏的概率为（1-（1-p1).^n(A_))
%  分解后，A有一定的概率会变成A_  (p2/timeInterval)  应该有p2>p3
%  也许，直接使用概率是错误的
beginA=50  ;   % 一开始有50个A
beginA_=0;     % 一开始没有坏的A
time=50000;    % 执行10000个时间段
p0=5;
p1=0.001;
p2=0.01;
p3=0.001;
result(1,1)=beginA;
result(1,2)=beginA_;
for i=1:time
    if i==40000
        p=1;
    end
    beginA=beginA+p0;
    [newA,newA_]=system2A(beginA,beginA_,p1,p2,p3);
    [newA_1]=system2A_(beginA,beginA_,p1);
    beginA=newA;
    beginA_=newA_+newA_1;
    result(i+1,1)=beginA;
    result(i+1,2)=beginA_;
end
hold on;
plot(result(:,1),'g');
plot(result(:,2),'r');
end
function [newA,newA_]=system2A(beginA,beginA_,p1,p2,p3)
newA=0;
newA_=0;
for i=1:beginA
    if rand>(1-(1-p1)^(beginA-1) )       
        if rand<p3+(1-(1-p1)^beginA_)*p2
            newA_=newA_+1;
        else
            newA=newA+1;
        end
    end
end
end
function [newA_1]=system2A_(beginA,beginA_,p1)
newA_1=0;
for i=1:beginA_
    if rand>1-(1-p1)^beginA
        newA_1=newA_1+1;
    end
end
end
        