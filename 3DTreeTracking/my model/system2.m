function result=system2()
%�������������һ�����ҷֽ��ϵͳA����������
%  ÿһ��ʱ�䣬��������һ������A ��p0 ��/timeInterval��
%  ͬʱA�Լ���һ���ĸ��ʻ�仵�����A_(ɥʧ����) ��p3 /timeInterval��
%  A���Էֽ⻵��A_  A+A_=A ,
%  ���������ʱ���ڣ��������������ĸ���Ϊp1(/timeInterval),��ô�������ĸ���Ϊ��1-��1-p1).^n(A_))
%  �ֽ��A��һ���ĸ��ʻ���A_  (p2/timeInterval)  Ӧ����p2>p3

% �任һ��˼·�����
% ÿһ��ʱ�䶼�����һ������A
% ��ЩA���ܻ��Լ��仵��p3
% ÿ��ʱ�䶼�����A+A_��ϸ���а���p1�ĸ���ѡ�� p1*C(A+A_)~ 2��ϸ���������
% ��������ΪA��������ɾ��1��������һ����A_�ĸ���Ϊp2
% ��������ΪA_
% һ��ΪA,��һ��ΪA_
beginA=5  ;   % һ��ʼ��50��A
beginA_=500;     % һ��ʼû�л���A
begin=[ones(beginA,1);zeros(beginA_,1)];
time=100000;    % ִ��10000��ʱ���
p0=5;
p1=0.01;  %ϸ�������ĸ���
p2=0.1;   %���ú�仵�ĸ���
p3=0.01;  %ֱ�ӱ仵�ĸ���
p4=0.5;    %�����������õĸ���
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
        