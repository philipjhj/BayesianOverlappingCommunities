function CompareTrueSampleZQ(obj,MAP,mymode)
% MAP = 0;
% mymode = 2;

zLabels = obj.zLabels;
qLabels = obj.qLabels;

if MAP
    [MAPP,MAPSampleZ,MAPSampleQ,iMax,jMax] = obj.findMAP;
    zMAP = obj.zSamples{MAPSampleZ};
    qMAP = obj.qSamples{MAPSampleQ};
    if mymode==1
        [zTrue, idx_o_Zt, idx_f_Zt] = sort_matrix(obj.zTrue,1);
        [qTrue, idx_o_Qt, idx_f_Qt] = sort_matrix(obj.qTrue,1);
    end
    [zMAP, idx_obj_Z, idx_fea_Z] = sort_matrix(zMAP,1);
    [qMAP, idx_obj_Q, idx_fea_Q] = sort_matrix(qMAP,1);
else
    zMAP=sum(cat(3,obj.zSamples{:}),3)/obj.T;
    qMAP=sum(cat(3,obj.qSamples{:}),3)/obj.T;
    if mymode == 1
        subplt = [223 224];
        subplot(2,2,1)
        imagesc(obj.zTrue)
        subplot(2,2,2)
        imagesc(obj.qTrue)
    else
        subplt = [211 212];
    end
    subplot(subplt(1))
    imagesc(zMAP,[0 1])
    set(gca,'YTick',1:size(zMAP,1));
    set(gca,'XTick',1:size(zMAP,2));
    if ~MAP
        idx_obj_Z = 1:size(zMAP,1);
        idx_fea_Z = 1:size(zMAP,2);
    end
    set(gca,'XTickLabel',idx_fea_Z);
    if ~isempty(zLabels)
        set(gca,'YTickLabel',zLabels(idx_obj_Z));
    else
        set(gca,'YTickLabel',idx_obj_Z);
    end
    
    %colorbar
    subplot(subplt(2))
    imagesc(qMAP,[0 1])
    set(gca,'YTick',1:size(qMAP,1));
    set(gca,'XTick',1:size(qMAP,2));
    if ~MAP
        idx_obj_Q = 1:size(qMAP,1);
        idx_fea_Q = 1:size(qMAP,2);
    end
    set(gca,'XTickLabel',idx_fea_Q);
    if ~isempty(qLabels)
        set(gca,'YTickLabel',qLabels(idx_obj_Q));
    else
        set(gca,'YTickLabel',idx_obj_Q);
    end
    %colorbar
    return
end


% Sorting to compare


% MIGHT BE USEFUL, DONT DELETE BEFORE REVISED (18/1/2016)
% colors=0;
% if colors == 1
%     last = find(sum(zTrue)<mean(sum(zTrue)),1)/2-1;
%     for i = 1:last
%         rows = zTrue(:,i)==1;
%         zTrue(rows,i) = zTrue(rows,i)+sum(zTrue(rows,i+1:last),2);
%         rows = qTrue(:,i)==1;
%         qTrue(rows,i) = qTrue(rows,i)+sum(qTrue(rows,i+1:last),2);
%
%         rows = zMAP(:,i)==1;
%         zMAP(rows,i) = zMAP(rows,i)+sum(zMAP(rows,i+1:last),2);
%         rows = qMAP(:,i)==1;
%         qMAP(rows,i) = qMAP(rows,i)+sum(qMAP(rows,i+1:last),2);
%     end
% end
% Plot all
% No of subplots to use
if mymode == 1
    subplt(1) = 221;
    subplt(3) = 222;
    subplt(2) = 223;
    subplt(4) = 224;
elseif mymode == 2
    subplt(2) = 121;
    subplt(4) = 122;
end


% Plot Z
if mymode == 1
    subplot(subplt(1))
    imagesc(zTrue)
    title('(a)')
    set(gca,'YTick',1:size(zTrue,1));
    set(gca,'XTick',1:size(zTrue,2));
    set(gca,'XTickLabel',idx_f_Zt);
    if ~isempty(zLabels)
        set(gca,'YTickLabel',zLabels(idx_o_Zt));
    else
        set(gca,'YTickLabel',idx_o_Zt);
    end
end
if any(mymode == [1 2])
    
    subplot(subplt(2))
    imagesc(zMAP,[0 1])
    set(gca,'YTick',1:size(zMAP,1));
    set(gca,'XTick',1:size(zMAP,2));
    set(gca,'XTickLabel',idx_fea_Z);
    if ~isempty(zLabels)
        set(gca,'YTickLabel',zLabels(idx_obj_Z));
    else
        set(gca,'YTickLabel',idx_obj_Z);
    end
    %set(gca,'YTickLabel',idx_obj_Z);
    %title(sprintf('MAP of Z at iteration: %d',MAPSampleZ))
    title('(c)')
end

% Plot Q
if mymode==1
    subplot(subplt(3))
    imagesc(qTrue)
    title('(b)')
    set(gca,'YTick',1:size(qTrue,1));
    set(gca,'XTick',1:size(qTrue,2));
    if ~isempty(qLabels)
        set(gca,'YTickLabel',qLabels);
    end
end
if any(mymode == [1 2])
    
    subplot(subplt(4))
    imagesc(qMAP,[0 1])
    %title(sprintf('MAP of Q at iteration: %d',MAPSampleQ))
    title('(d)')
    set(gca,'YTick',1:size(qMAP,1));
    set(gca,'XTick',1:size(qMAP,2));
    set(gca,'XTickLabel',idx_fea_Q);
    if ~isempty(qLabels)
        set(gca,'YTickLabel',qLabels(idx_obj_Q));
    else
        set(gca,'YTickLabel',idx_obj_Q);
    end
    %set(gca,'YTickLabel',idx_obj_Q);
end
end

% function [M, idx_obj, idx_fea] = sort_matrix(M,iter,idx_in)
%
% if iter == 1
%     sumMc=sum(M,1);
%     [~,idxMc]=sort(sumMc,'descend');
%     idx_fea = idxMc;
%     M = M(:,idxMc);
% end
%
% % Order of sorting
% if any(iter==[1 3])
%     sort_string = 'descend';
% elseif any(iter ==[2])
%     sort_string = 'ascend';
% end
%
% % Sort rows of current matrix
% [~,idxMr]=sort(M(:,1),sort_string);
% M = M(idxMr,:);
%
% % If first round
% if iter == 1
%     idx_in = (1:size(M,1))';%idxMr;
% end
%
% % Assign intial indices to current sub-matrix
% idx_temp = idx_in(idxMr);
%
% % Only sort if there is more than one row
% if size(M,2) > 1 && size(M,1) > 1
%     % To account for order of sorting
%     if any(iter == [ 1 3 ])
%         upper_m = M(:,1)>0;
%         lower_m = (M(:,1)==0);
%     else
%         upper_m = M(:,1)==0;
%         lower_m = (M(:,1)>0);
%     end
%
%     idxMr_up = idx_temp(upper_m);
%     idxMr_lo = idx_temp(lower_m);
%
%     M_1 = M(upper_m,2:end);
%     [M_1, idx_1] = sort_matrix(M_1,2,idxMr_up);
%     M_2 = M(lower_m,2:end);
%     [M_2, idx_2] = sort_matrix(M_2,3,idxMr_lo);
%     M12 = [M_1 ; M_2];
%     M = [M(:,1) M12];
%     idx_obj = [idx_1; idx_2];
% else
%     idx_obj = idx_in(idxMr);
% end
%
% end