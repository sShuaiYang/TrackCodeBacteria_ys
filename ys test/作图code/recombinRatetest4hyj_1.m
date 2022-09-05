function [out] = recombinRatetest4hyj_1(A1,A2)
%������������ΪA1������ΪA2����׼��B1��B2Ϊ���������1%��50%�仯��������̬�ֲ�
%�����������
%���������Բ���,
m=3;%Ϊʵ���ظ�����
N=200;%�ʱ�׼��ı仯�ĸ�����Ĭ��100
B1=linspace(A1*0.01,A1*0.5,N);
B2=linspace(A2*0.01,A2*0.5,N);
k=1;
for i=1:N
    for j=1:N
        
        seq_p=[];
        for h=1:100
            %�ظ�һ�ٴμ���ƽ��pֵ
            seqA1=random('Normal',A1,B1(i),m,1);
            seqA2=random('Normal',A2,B2(j),m,1);
            [~,seq_p(h),~,~] = ttest2(seqA1,seqA2);
        end
        out(1,k)=B1(i);
        out(2,k)=B2(j);
        out(3,k)=mean(seq_p);
        k=k+1;
    end
    
    
end

a1=out(3,:)<0.001;
out1=out(:,a1);
a2=out(3,:)>=0.001&out(3,:)<0.05;
out2=out(:,a2);
a3=out(3,:)>=0.05;
out3=out(:,a3);

figure, plot3(out1(1,:)/A1,out1(2,:)/A2,out1(3,:),'Marker','.','LineStyle','none','Color',[44,127,184]/255)
hold on, plot3(out2(1,:)/A1,out2(2,:)/A2,out2(3,:),'Marker','.','LineStyle','none','Color',[127,205,187]/255)
plot3(out3(1,:)/A1,out3(2,:)/A2,out3(3,:),'Marker','.','LineStyle','none','Color',[237,248,177]/255)
% ���� zlabel
zlabel('ttest-p');

% ���� ylabel
ylabel('sigma%A2');

% ���� xlabel
xlabel('sigma%A1');

[X,Y] = meshgrid(B1,B2);
Z = reshape(out(3,:),[N,N]);
figure,surf(X/5,Y/10,Z,'EdgeColor','none')
view(2)
colormap(jet)


end
