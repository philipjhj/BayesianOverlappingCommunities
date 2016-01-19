%function logPversusHistogramofSamples
col = hsv(8);

for i = 1:4
    for j = 1:2
        mySuffStats = SuffStatsClass;
        current_pms = myResults.pmsSamples{end};
        logMid=log(current_pms(i,j))/log(10);
        test_values=logspace(logMid-0.1, logMid+0.1);
        pmsAll = repmat(current_pms,1,1,numel(test_values));
        pmsAll(i,j,:) = test_values;
        
        for k = 1:numel(test_values)
            logPs(k)=LogJoint(true,A,myResults.zSamples{end},myResults.qSamples{end},pmsAll(:,:,k),mySuffStats);
        end
        
        %logPfun = @(x) exp(IntegralFunctionLogP(x,A,zTrue,qTrue,pms,mySuffStats,i,j))
        %area = integral(logPfun,0,10)
        
        
        subplot(2,4,(j-1)*4+i);
        yy = exp(logPs)./(trapz(exp(logPs)));
        yy = logPs;
        plot((test_values),yy,'linewidth',2,'color',col(sub2ind([2,4],j,i),:))
        p=0.5;
        %xlim([logMid*(1-p) (logMid)*(1+p)])
        test_values;
        %break;
%         hold on
%         pmsSamples=cat(3,myResults.pmsSamples{966:end});
%         %trapz(exp(logP))
%         %sum(exp(logP))
%         
%         myvalues = reshape(pmsSamples(i,j,:),1,numel(pmsSamples(i,j,:)));
%         histogram(myvalues,[0:xres:xmax],'Normalization','pdf')
%         hold off
    end
    %break;
end