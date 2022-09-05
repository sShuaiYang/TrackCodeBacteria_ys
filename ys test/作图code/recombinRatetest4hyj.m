function [a,b,A,B,p_data] = recombinRatetest4hyj(A1,B1,A2,B2,N)
%期望为A1，标准差为B1，期望为A2， 标准差为B2，生成正态分布的两组随机数
%N 为数目 计算显著性差异
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

