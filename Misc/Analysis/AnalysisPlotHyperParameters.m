pms_data = cat(3,myResults.pmsSamples{:});

col = hsv(8);
range = 1:numel(pms_data(1,1,:));
for i = 1:size(myResults.pms,1)
    for j = 1:size(myResults.pms,2)
        plot(log(reshape(pms_data(i,j,range),1,numel(pms_data(i,j,range)))),'color',col(sub2ind([2,4],j,i),:),'linewidth',1)
        hold on
    end
end
legend on