%function logPversusHistogramofSamples
col = hsv(16);
logpcontri = zeros(4,2,numel(test_values));
for i = 1:4
    for j = 1:2
        mySuffStats = SuffStatsClass;
        
        current_pms = myResults.pmsSamples{end};
        logMid=log(current_pms(i,j))/log(10);
        test_values=logspace(logMid-0.1, logMid+0.1);
        %test_values_log = logMid+[-logspace(0,0.01) +logspace(0,0.01)];
        %test_values = exp(test_values_log);
        %i=4;
        %j=2;
        pmsAll = repmat(current_pms,1,1,numel(test_values));
        pmsAll(i,j,:) = test_values;
        
        for k = 1:numel(test_values)
            [~,~,newsuff]=LogJoint(true,A,myResults.zSamples{end},myResults.qSamples{end},pmsAll(:,:,k),mySuffStats);
            
            if i == 1
                %%
                %logpcontri(i,1,k) = betaln(newsuff.NZp,newsuff.NZm) - betaln(pmsAll(1,1,k),pmsAll(1,2,k));
                
                %%
                np = newsuff.NZp - pmsAll(1,1,1);
                nm = newsuff.NZm - pmsAll(1,2,1);
                bp = permute( pmsAll(1,1,:),[1,3,2]);
                bm = zeros(size(bp))+pmsAll(1,2,1); 
                
                np = np(end:end);
                nm = nm(end:end);
                
                bmax = 100; 
                K = 100; 
                x = linspace(0.01,bmax, K); 
                
                [BP,BM] = meshgrid(x,x,K);
                close all;
                %%
                figure
                %%
                bb = linspace(1,100,100);
                plot(bb, betaln(2+bb,4
+bb)-betaln(bb,bb));
                %%
                y = sum(bsxfun(@minus,betaln(bsxfun(@plus, np,BP(:)'), bsxfun(@plus,nm,BM(:)')),betaln(BP(:)',BM(:)')),1);
                y = reshape(y,[K,K]);
                [~,ii] = max(y(:));
                [I,J] = ind2sub(size(y), ii);
                
                surf(BP,BM,y);
                [v,i] = max(y);
                [v,j] = max(v);
                [i(j),j, I, J]
                
                plot(max(y))
                
                
                
                %%
                y = sum(bsxfun(@minus,betaln(bsxfun(@plus, np,bp), bsxfun(@plus,nm,bm)),betaln(bp,bm)),1);
                y = squeeze(y);
                semilogx(bp, y);
                %%
                
                betaln(newsuff.NZp,newsuff.NZm) - betaln(pmsAll(1,1,k),pmsAll(1,2,k));
                
%                 logpcontri(i,2,k) = size(myResults.zTrue,2)*;
            elseif i == 2
                logpcontri(i,1,k) = sum(betaln(newsuff.NQp,newsuff.NQm));
                logpcontri(i,2,k) = size(myResults.qTrue,2)*betaln(pmsAll(2,1,k),pmsAll(2,2,k));
            elseif i == 3
                logpcontri(i,2,k) = size(A,2)*betaln(pmsAll(3,1,k),pmsAll(3,2,k));
                %logpcontri(i,2,k) = betaln(newsuff.NPaps,newsuff.NPams);
                logpcontri(i,1,k) = sum(betaintegral_noverbose(newsuff.NPaps,newsuff.NPams,newsuff.NPbps,newsuff.NPbms,50000));
            else
                logpcontri(i,2,k) = size(A,2)*betaln(pmsAll(4,1,k),pmsAll(4,2,k));
                %logpcontri(i,2,k) = betaln(newsuff.NPbps,newsuff.NPbms);
                logpcontri(i,1,k) = sum(betaintegral_noverbose(newsuff.NPaps,newsuff.NPams,newsuff.NPbps,newsuff.NPbms,50000));
            end
        end
        
        subplot(2,4,(j-1)*4+i);
        for m = 1:2
            plot(log(test_values),reshape(logpcontri(i,m,:),1,numel(test_values)),'linewidth',2,'color',col(sub2ind([2,2,4],m,j,i),:))
            hold on
        end
        hold off
    end
    %break;
end