figure
col=hsv(16)
a = ones(16,100)
b=1:16

c=repmat(b',1,100).*a;

for i = 1:16
    plot(c(i,:),'color',col(i,:),'linewidth',3)
    hold on
end