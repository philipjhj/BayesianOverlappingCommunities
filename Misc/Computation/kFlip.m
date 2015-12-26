function [rand_idx,permutations,perm_changes,full_shadow,perm_order]=k_flip2(matrix1,matrix2_row,k)
% Given a matrix1, k random indices are sampled and
% the order of "easiest-to-calculate" permutations is found
% for these indices in the opposit matrix (matrix2_row)
% I.e if Q is given, permutations of the random indices are found for Z,
% and vice versa

if k == 1
    error('k must be 2 or more')
end

d=size(matrix1,2);
rand_idx = sort(randsample(d,k));
%rand_idx=[1 3 4];

not_selected=[];
for i = 1:length(matrix2_row)
    if all(i~=rand_idx) && matrix2_row(i) == 1
        not_selected=[not_selected i];
    end
end

matrix2_shadow = sum(matrix1(:,find(matrix2_row)),2)>0;

init_shadow=sum(matrix1(:,rand_idx),2)>0;

[~,shadow_first] = find_diff_order((matrix2_shadow)',init_shadow',rand_idx);

[perm_order,full_shadow]=find_perm_order(matrix1,matrix2_row,rand_idx,k,{},not_selected,[]);

full_shadow = [shadow_first full_shadow];

permutations = zeros(2^k,k);
permutations(1,:) = ones(1,k);

reset_points=[1];
for i = 1:k-1
    reset_points = [reset_points 1+reset_points(2^(i-1))*2  reset_points];
end

 perm_order;
for i = 1:length(perm_order)
    permutations(i+1,:) = permutations(i,:);
    permutations(i+1,find(rand_idx==perm_order(i))) = 0;
    if mod(i,2) == 0
        permutations(i+1,ismember(rand_idx,unique(perm_order((i-reset_points(i/2)):(i-1))))) = 1;
    end
end

perm_changes = cell(1,2^k);
perm_changes{1} = rand_idx(find(abs(permutations(1,:)-matrix2_row(rand_idx))));
for i = 1:length(permutations)-1
    perm_changes{i+1} = rand_idx(find(abs(permutations(i,:)-permutations(i+1,:))));
end

for i = 1:length(full_shadow)
    full_shadow{i} = unique(full_shadow{i});
end


end

% Finds the complete list of which indices to flip in order
function [perm_order,full_shadow,rand_idx] = find_perm_order(matrix1,matrix2_row,rand_idx,k,full_shadow,not_selected,perm_order_prev);

matrix2_shadow = sum(matrix1(:,not_selected),2)>0;

[diff_order,full_shadow_test] = find_diff_order(matrix2_shadow',matrix1(:,rand_idx)',rand_idx);
%[~,full_shadow_test] = find_diff_order(matrix2_shadow',matrix1(:,rand_idx)',rand_idx);

[perm_order,full_shadow] = add_changes(diff_order(1:end),matrix1,[],full_shadow,full_shadow_test,matrix2_shadow);

if length(diff_order) > 2
    perm_order = [perm_order diff_order(1)];
    
    [~,next_shadow] = find_diff_order((matrix2_shadow+sum(matrix1(:,diff_order(2:end)),2))'>0,matrix1(:,diff_order(1))',diff_order(1));
    full_shadow{end+1} = [next_shadow{1} full_shadow_test{2:end}];
    
    perm_order_prev = [perm_order_prev perm_order];
    
    [new_perms,full_shadow] = find_perm_order(matrix1,matrix2_row,diff_order(2:end),k,full_shadow,not_selected,perm_order_prev);
    perm_order = [perm_order new_perms];
    
else
    [~,next_shadow] = find_diff_order(matrix2_shadow',matrix1(:,diff_order(2))',diff_order(2));
    full_shadow{end+1} =  [full_shadow_test{1} next_shadow{end}];
    full_shadow{end+1} = next_shadow{end};
    
    perm_order = [perm_order diff_order];
end
end

% Function adds which indices to flip in order given current position in
% "tree"
function [perm_order,full_shadow]=add_changes(diff_order,matrix,perm_order,full_shadow,full_shadow_test,matrix2_shadow)

% cellfun(@length,full_shadow_test);

n_full_covered=sum(cellfun(@length,full_shadow_test)==0);
idx_full_covered=find(cellfun(@length,full_shadow_test)==0);

if n_full_covered > 2
%     length(full_shadow_test);
    full_covered=[1];
    for i = 1:n_full_covered-1
        full_covered = [full_covered i+1 full_covered];
    end
%     full_covered;
%     idx_full_covered(full_covered);
%     diff_order;
%     diff_order(idx_full_covered(full_covered));
    perm_order = [perm_order diff_order(idx_full_covered(full_covered))];
    
    for i = 1:length(full_covered)
        full_shadow{(end+1)} = [];
    end
    
else
    n_full_covered=0;
end


for i = (1+n_full_covered):length(diff_order)-1
    if i < 3
        perm_order=[perm_order diff_order((end-i+1):end)];
        if i == 1
            full_shadow{end+1}=full_shadow_test{end};
        elseif i == 2
            
            test_shadow = matrix(:,[diff_order(1:(end-i))])>0;
            if length(diff_order(2:end))>2
                test_shadow=sum(test_shadow,2);
            end
            test_shadow=(test_shadow+matrix2_shadow)>0;
            
            [~,next_shadow] = find_diff_order(test_shadow',matrix(:,diff_order(end))',diff_order(end));
            full_shadow{end+1} = [full_shadow_test{end-i+1} next_shadow{end}];
            full_shadow{end+1} = next_shadow{end};
            
        end
    else
        perm_order = [perm_order diff_order(end+1-i)];
        
        new_shadow = (matrix(:,diff_order(1))+matrix2_shadow)>0;
        tmp_diff_order=diff_order(2:end);
        tmp_diff_order(end+1-i) = [];
        [new_diff_order,new_test_shadow]=find_diff_order(new_shadow',matrix(:,tmp_diff_order)',tmp_diff_order);
        
        more_shadow=cat(2,new_test_shadow{(end-i+2):end});
        
        
        full_shadow{end+1} = [full_shadow_test{end-i+1} more_shadow];
%         i
%         length([diff_order(1) new_diff_order((end-i+2):end)])
%         length([full_shadow_test{1} new_test_shadow((end-i+2):end)])
%         [full_shadow_test{1} new_test_shadow((end-i+2):end)]
        [perm_order,full_shadow]=add_changes([diff_order(1) new_diff_order((end-i+2):end)],matrix,perm_order,...
            full_shadow,[full_shadow_test{1} new_test_shadow((end-i+2):end)],matrix2_shadow);
    end
end
end

% diff_order is in order of most different from highest level active node
function [diff_order,test_shadow] = find_diff_order(shadow,remaining,test_idx,test_shadow)
reps=size(remaining,1);
most_shadow = sum((repmat(shadow,reps,1)-remaining)==-1,2);
[values,order]=sort(most_shadow,'descend');

diff_order=test_idx(order(1));

if ~exist('test_shadow','var')
    test_shadow = {};
end

if values(1)==0
    test_shadow{end+1}=[];
else
    shadow_diff=shadow-remaining(order(1),:);
    test_shadow{end+1}=find(shadow_diff<0);
end

if size(remaining,1)>1
    new_shadow = (shadow+remaining(order(1),:))>0;
    remaining(order(1),:) = [];
    test_idx(order(1)) = [];
    [add_order,test_shadow]=find_diff_order(new_shadow,remaining,test_idx,test_shadow);
    diff_order = [diff_order add_order];
end
end
