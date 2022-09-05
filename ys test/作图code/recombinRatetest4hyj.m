function [a,b,A,B,p_data] = recombinRatetest4hyj(A1,B1,A2,B2,N)
%����ΪA1����׼��ΪB1������ΪA2�� ��׼��ΪB2��������̬�ֲ������������
%N Ϊ��Ŀ ���������Բ���
a = random('Normal',A1,B1,N,1);
b = random('Normal',A2,B2,N,1);
A=nchoosek(a,3);
B=nchoosek(b,3);
p_data=[];
for i=1:size(A,1)
    for j=1:size(B,1)
        [~,p_data(i,j),~,~] = ttest2(A(i,:)',B(i,:)');
    end
end
% figure, histogram(a)
% hold on,histogram(b)
% 
% figure, histogram(p)

end

