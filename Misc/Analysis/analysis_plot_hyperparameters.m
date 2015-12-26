pms_data = cat(3,pms_all{:});


col = hsv(8);
range = 1200:numel(pms_data(1,1,:));
for i = 1:size(pms,1)
    for j = 1:size(pms,2)
        plot((reshape(pms_data(i,j,range),1,numel(pms_data(i,j,range)))),'color',col(sub2ind([4,2],i,j),:))
        hold on
    end
end