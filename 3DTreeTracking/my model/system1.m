function result=system1()
%�������������һ�����ҷֽ��ϵͳA����������
%  ÿһ��ʱ�䣬��������һ������A ��p0 ��/timeInterval��
%  ͬʱA�Լ���һ���ĸ��ʻ�仵�����A_(ɥʧ����) ��p3 /timeInterval��
%  A���Էֽ⻵��A_  A+A_=A ,
%  ���������ʱ���ڣ��������������ĸ���Ϊp1(/timeInterval),��ô�������ĸ���Ϊ��1-��1-p1).^n(A_))
%  �ֽ��A��һ���ĸ��ʻ���A_  (p2/timeInterval)  Ӧ����p2>p3
%  Ҳ��ֱ��ʹ�ø����Ǵ����
beginA=50  ;   % һ��ʼ��50��A
beginA_=0;     % һ��ʼû�л���A
time=50000;    % ִ��10000��ʱ���
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
        