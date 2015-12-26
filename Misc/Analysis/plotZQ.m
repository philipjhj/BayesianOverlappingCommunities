function [idx_obj_Z, idx_obj_Q,idx_fea_Z,idx_fea_Q] = plotZQ(Z,Q,colors,mode,Z_labels,Q_labels)
% Visualize the clusters Z and Q
% Orginally and Sorted
%   Colors: 
%   Mode:
%               1: Plot all matrices
%               2: Plot only Z
%               3: Plot only W

Z_org = Z;
Q_org = Q;

[Z, idx_obj_Z, idx_fea_Z] = sort_matrix(Z,1);
[Q, idx_obj_Q, idx_fea_Q] = sort_matrix(Q,1);

% If sorting is messing something up
if 0 ~= sum(sort(sum(Z_org,1),'descend')-sum(Z,1))+sum(sort(sum(Q_org,1),'descend')-sum(Q,1))
    display('Sorting is wrong');
end


%
if colors == 1
    last = find(sum(Z)<mean(sum(Z)),1)/2-1;
    for i = 1:last
        rows = Z(:,i)==1;
        Z(rows,i) = Z(rows,i)+sum(Z(rows,i+1:last),2);
    end
    
    for i = 1:last %size(Q,2)
        rows = Q(:,i)==1;
        Q(rows,i) = Q(rows,i)+sum(Q(rows,i+1:last),2);
    end
end

% No of subplots to use
if mode == 1
    subplt(1) = 221;
    subplt(3) = 222;
    subplt(2) = 223;
    subplt(4) = 224;
elseif mode == 2
    subplt(1) = 211;
    subplt(2) = 212;
elseif mode == 3
    subplt(3) = 211;
    subplt(4) = 212;
end

% Plot Z
if any(mode == [1 2])
subplot(subplt(1))
imagesc(Z_org)
title('Z')
set(gca,'YTick',1:size(Z,1));
set(gca,'XTick',1:size(Z,2));
if exist('Z_labels')
    set(gca,'YTickLabel',Z_labels);
end



subplot(subplt(2))
imagesc(Z)
set(gca,'YTick',1:size(Z,1));
set(gca,'XTick',1:size(Z,2));
set(gca,'XTickLabel',idx_fea_Z);
if exist('Z_labels')
    set(gca,'YTickLabel',Z_labels(idx_obj_Z));
end
%set(gca,'YTickLabel',idx_obj_Z);
title('Z - Sorted')
end

% Plot Q
if any(mode==[1 3])
subplot(subplt(3))
imagesc(Q_org)
title('Q')
set(gca,'YTick',1:size(Q,1));
set(gca,'XTick',1:size(Q,2));
if exist('Q_labels')
    set(gca,'YTickLabel',Q_labels);
end


subplot(subplt(4))
imagesc(Q)
title('Q - Sorted')
set(gca,'YTick',1:size(Q,1));
set(gca,'XTick',1:size(Q,2));
set(gca,'XTickLabel',idx_fea_Q);
if exist('Q_labels')
    set(gca,'YTickLabel',Q_labels(idx_obj_Q));
end
%set(gca,'YTickLabel',idx_obj_Q);
end
end

