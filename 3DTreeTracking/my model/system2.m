function result=system2()
%　定义的是这样一个自我分解的系统A，规则如下
%  每一段时间，会生产出一定量的A （p0 个/timeInterval）
%  同时A自己有一定的概率会变坏，变成A_(丧失功能) （p3 /timeInterval）
%  A可以分解坏的A_  A+A_=A ,
%  假设在这段时间内，两个蛋白相遇的概率为p1(/timeInterval),那么遇到坏的概率为（1-（1-p1).^n(A_))
%  分解后，A有一定的概率会变成A_  (p2/timeInterval)  应该有p2>p3

% 变换一下思路，变成
% 每一段时间都会产生一定量的A
% 这些A可能会自己变坏，p3
% 每段时间都会从这A+A_个细菌中按照p1的概率选出 p1*C(A+A_)~ 2组细菌进行配对
% 若两个都为A，则随意删除1个，另外一个变A_的概率为p2
% 若两个都为A_
% 一个为A,另一个为A_
beginA=5  ;   % 一开始有50个A
beginA_=500;     % 一开始没有坏的A
begin=[ones(beginA,1);zeros(beginA_,1)];
time=100000;    % 执行10000个时间段
p0=5;
p1=0.01;  %细菌相遇的概率
p2=0.1;   %作用后变坏的概率
p3=0.01;  %直接变坏的概率
p4=0.5;    %相遇后发生作用的概率
result(1,1)=beginA;
result(1,2)=beginA_;
for i=1:time
    disp(i)
    if rand<=0.5
        begin=[begin;ones(fix(p0*numel(begin(begin==1))/numel(begin)),1)];
    else
        begin=[begin;ones(ceil(p0*numel(begin(begin==1))/numel(begin)),1)];
    end
    begin=systemAapoptosis(begin,p3);
    begin=system2Ameet(begin,p1,p2,p4);
    result(i+1,1)=numel(begin(begin==1));
    result(i+1,2)=numel(begin(begin==0));
end
hold on;
plot(result(:,1),'g');
plot(result(:,2),'r');
end
function begin=systemAapoptosis(begin,p3)
for i=1:numel(begin)
    if begin(i)==1 && rand<p3
        begin(i)=0;
    end
end
end
function [begin]=system2Ameet(begin,p1,p2,p4)
meetPair=numel(begin)*(numel(begin)-1)/2*p1;
for i=1:meetPair
    num1=ceil(numel(begin)*rand);
    num2=ceil(numel(begin)*rand);
    while num2==num1
        num2=ceil(numel(begin)*rand);
    end
    if begin(num1)==0 && begin(num2)==0
        continue
    end
    if begin(num1)==1 && begin(num2)==0
        if rand<p4
            if rand<p2
                begin(num1)=0;
            end
            begin(num2)=[];
        end
        continue
    end
    if begin(num1)==0 && begin(num2)==1
        if rand<p4
            if rand<p2
                begin(num2)=0;
            end
            begin(num1)=[];
        end
        continue
    end
    if begin(num1)==1 && begin(num2)==1
        if rand<p4
            if rand<=0.5
                if rand<p2
                    begin(num1)=0;
                end
                begin(num2)=[];
            else
                if rand<p2
                    begin(num2)=0;
                end
                begin(num1)=[];
            end
        end
        continue
    end
end
end
        