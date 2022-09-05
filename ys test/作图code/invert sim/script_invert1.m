figure,
color=jet(21);
for i=1:21
    subplot(2,3,1)
    plot(simdata{i}.t,simdata{i}.x(:,2),'color',color(i,:))
    hold on    
end
xlabel('t (s)')
ylabel('mRNA1')

for i=1:21
    subplot(2,3,2)
    plot(simdata{i}.t,simdata{i}.x(:,3),'color',color(i,:))
    hold on    
end
xlabel('t (s)')
ylabel('protein1')

for i=1:21
    subplot(2,3,3)
    plot(simdata{i}.t,simdata{i}.x(:,4),'color',color(i,:))
    hold on    
end
xlabel('t (s)')
ylabel('DNA2')

for i=1:21
    subplot(2,3,4)
    plot(simdata{i}.t,simdata{i}.x(:,5),'color',color(i,:))
    hold on    
end
xlabel('t (s)')
ylabel('mRNA2')

for i=1:21
    subplot(2,3,5)
    plot(simdata{i}.t,simdata{i}.x(:,6),'color',color(i,:))
    hold on    
end
xlabel('t (s)')
ylabel('protein2')

for i=1:21
    subplot(2,3,6)
    plot(simdata{i}.t,simdata{i}.x(:,7),'color',color(i,:))
    hold on    
end
xlabel('t (s)')
ylabel('DNA2-protein1')

for j=1:21
    p1(j)=simdata{j}.x(end,3);    
    p2(j)=simdata{j}.x(end,6);
end

for j=1:21
    x1(j)=simdata{j}.x(end,2);    
    
end


