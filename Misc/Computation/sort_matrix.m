function [M, idx_obj, idx_fea] = sort_matrix(M,iter,idx_in)

if iter == 1
    sumMc=sum(M,1);
    [~,idxMc]=sort(sumMc,'descend');
    idx_fea = idxMc;
    M = M(:,idxMc);
end

% Order of sorting
if any(iter==[1 3])
    sort_string = 'descend';
elseif any(iter ==[2])
    sort_string = 'ascend';
end

% Sort rows of current matrix
[~,idxMr]=sort(M(:,1),sort_string);
M = M(idxMr,:);

% If first round
if iter == 1
    idx_in = (1:size(M,1))';%idxMr;
end

% Assign intial indices to current sub-matrix
idx_temp = idx_in(idxMr);

% Only sort if there is more than one row
if size(M,2) > 1 && size(M,1) > 1
    % To account for order of sorting
    if any(iter == [ 1 3 ])
        upper_m = M(:,1)>0;
        lower_m = (M(:,1)==0);
    else
        upper_m = M(:,1)==0;
        lower_m = (M(:,1)>0);
    end
    
    idxMr_up = idx_temp(upper_m);
    idxMr_lo = idx_temp(lower_m);
    
    M_1 = M(upper_m,2:end);
    [M_1, idx_1] = sort_matrix(M_1,2,idxMr_up);
    M_2 = M(lower_m,2:end);
    [M_2, idx_2] = sort_matrix(M_2,3,idxMr_lo);
    M12 = [M_1 ; M_2];
    M = [M(:,1) M12];
    idx_obj = [idx_1; idx_2];
else
    idx_obj = idx_in(idxMr);
end

end