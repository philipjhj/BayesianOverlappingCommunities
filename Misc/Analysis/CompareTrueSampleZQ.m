function CompareTrueSampleZQ(zTrue,qTrue,zSample,qSample)

% Saving Original Ordering
% zTrueOrg = zTrue;
% qTrueOrg = qTrue;
% zSampleOrg = zSample;
% qSampleOrg = qSample;

% Sorting to compare
[zTrue, idx_obj_Z, idx_fea_Z] = sort_matrix(zTrue,1);
[qTrue, idx_obj_Q, idx_fea_Q] = sort_matrix(qTrue,1);
[zSample, idx_obj_Z, idx_fea_Z] = sort_matrix(zSample,1);
[qSample, idx_obj_Q, idx_fea_Q] = sort_matrix(qSample,1);

% If sorting is messing something up
% if 0 ~= sum(sort(sum(zTrueOrg,1),'descend')-sum(zTrue,1))+sum(sort(sum(qTrueOrg,1),'descend')-sum(qTrue,1))
%     display('Sorting is wrong');
% end
% if 0 ~= sum(sort(sum(zSampleOrg,1),'descend')-sum(zSample,1))+sum(sort(sum(qSampleOrg,1),'descend')-sum(qTrue,1))
%     display('Sorting is wrong');
% end

colors=0;
if colors == 1
    last = find(sum(zTrue)<mean(sum(zTrue)),1)/2-1;
    for i = 1:last
        rows = zTrue(:,i)==1;
        zTrue(rows,i) = zTrue(rows,i)+sum(zTrue(rows,i+1:last),2);
        rows = qTrue(:,i)==1;
        qTrue(rows,i) = qTrue(rows,i)+sum(qTrue(rows,i+1:last),2);
        
        rows = zSample(:,i)==1;
        zSample(rows,i) = zSample(rows,i)+sum(zSample(rows,i+1:last),2);
        rows = qSample(:,i)==1;
        qSample(rows,i) = qSample(rows,i)+sum(qSample(rows,i+1:last),2);
    end
end

mymode = 1; % Plot all
% No of subplots to use
if mymode == 1
    subplt(1) = 221;
    subplt(3) = 222;
    subplt(2) = 223;
    subplt(4) = 224;
elseif mymode == 2
    subplt(1) = 211;
    subplt(2) = 212;
elseif mymode == 3
    subplt(3) = 211;
    subplt(4) = 212;
end


% Plot Z
if any(mymode == [1 2])
subplot(subplt(1))
imagesc(zTrue)
title('True Z')
set(gca,'YTick',1:size(zTrue,1));
set(gca,'XTick',1:size(zTrue,2));
if exist('Z_labels')
    set(gca,'YTickLabel',Z_labels);
end

subplot(subplt(2))
imagesc(zSample)
set(gca,'YTick',1:size(zSample,1));
set(gca,'XTick',1:size(zSample,2));
set(gca,'XTickLabel',idx_fea_Z);
if exist('Z_labels')
    set(gca,'YTickLabel',Z_labels(idx_obj_Z));
end
%set(gca,'YTickLabel',idx_obj_Z);
title('Sampled Z')
end

% Plot Q
if any(mymode==[1 3])
subplot(subplt(3))
imagesc(qTrue)
title('True Q')
set(gca,'YTick',1:size(qTrue,1));
set(gca,'XTick',1:size(qTrue,2));
if exist('Q_labels')
    set(gca,'YTickLabel',Q_labels);
end


subplot(subplt(4))
imagesc(qSample)
title('Sampled Q')
set(gca,'YTick',1:size(qTrue,1));
set(gca,'XTick',1:size(qTrue,2));
set(gca,'XTickLabel',idx_fea_Q);
if exist('Q_labels')
    set(gca,'YTickLabel',Q_labels(idx_obj_Q));
end
%set(gca,'YTickLabel',idx_obj_Q);
end


end