function [rand_idx,permutations,perm_order,full_shadow]=k_flip(matrix1,matrix2_row,k)
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

not_selected=[];
for i = 1:length(matrix2_row)
    if all(i~=rand_idx) && matrix2_row(i) == 1
        not_selected=[not_selected i];
    end
end


[perm_order,full_shadow]=find_perm_order(matrix1,matrix2_row,rand_idx,k,{},not_selected,[]);

% test diff_order
%find_diff_order([1 1 0 0 0], [0 0 1 0 0; 1 1 0 0 0; 0 0 1 1 0],2:4)

permutations(1,:) = ones(1,k);
for i = 1:length(perm_order)
    permutations(i+1,:) = permutations(i,:);
    permutations(i+1,find(rand_idx==perm_order(i)))=abs(permutations(i+1,find(rand_idx==perm_order(i)))-1);
end
end

% Finds the complete list of which indices to flip in order
function [perm_order,full_shadow,rand_idx] = find_perm_order(matrix1,matrix2_row,rand_idx,k,full_shadow,not_selected,perm_order_prev);
% most_active = sum(matrix1(:,rand_idx),1);
% [~, order] = sort(most_active,'descend');
% biggest_active = rand_idx(order(1));
% remaining = rand_idx(order(2:end));

%rand_idx
%not_selected

if length(perm_order_prev) > 0
    actives=find_actives(perm_order_prev);
end


matrix2_shadow = sum(matrix1(:,not_selected),2)>0;

%[diff_order,full_shadow_test] = find_diff_order(matrix1(:,biggest_active)',matrix1(:,remaining)',remaining);
[diff_order,~] = find_diff_order(matrix2_shadow',matrix1(:,rand_idx)',rand_idx);
[~,full_shadow_test] = find_diff_order(matrix2_shadow',matrix1(:,rand_idx)',rand_idx);

[perm_order,full_shadow] = add_changes(diff_order(2:end),diff_order(1),matrix1,[],full_shadow,full_shadow_test,matrix2_shadow);
%[perm_order,full_shadow] = add_changes(remaining,biggest_active,matrix1,[],full_shadow,full_shadow_test);


if length(diff_order) > 2
    full_shadow{end+1} = full_shadow_test{1};
    perm_order = [perm_order diff_order(1)];
    perm_order_prev = [perm_order_prev perm_order];
    [new_perms,full_shadow] = find_perm_order(matrix1,matrix2_row,diff_order(2:end),k,full_shadow,not_selected,perm_order_prev);
    %[new_perms,full_shadow] = find_perm_order(matrix1,remaining,k,full_shadow);
    %full_shadow{end+1: = [full_shadow 0 more_shadow];
    perm_order = [perm_order new_perms];
    %perm_order = [perm_order biggest_active new_perms];
else
    full_shadow{end+1} =  full_shadow_test{1};
    %
    %     test_shadow = matrix(:,[ biggest_active diff_order(1:(end-i))])>0
    %     if length(diff_order)>2
    %         test_shadow=sum(test_shadow)
    %     end
    %     %test_shadow = sum(matrix(:,[ biggest_active diff_order(1:(end-i))]))'>0
    %     [~,next_shadow] = find_diff_order(test_shadow',matrix(:,diff_order(end))',diff_order(end));
    %     full_shadow{end+1} = next_shadow{end};
    %
    full_shadow{end+1} = full_shadow_test{2};%full_shadow_test{end};
    perm_order = [perm_order diff_order];
    %perm_order = [perm_order biggest_active remaining];
end
end

% Function adds which indices to flip in order given current position in
% "tree"
function [perm_order,full_shadow]=add_changes(diff_order,biggest_active,matrix,perm_order,full_shadow,full_shadow_test,matrix2_shadow)
for i = 1:length(diff_order)
    if i < 3
        perm_order=[perm_order diff_order((end-i+1):end)];
        if i < 2
            full_shadow{end+1}=full_shadow_test{end};
        else
            full_shadow{end+1} = full_shadow_test{end-i+1};
            test_shadow = matrix(:,[ biggest_active diff_order(1:(end-i))])>0;
            if length(diff_order)>2
                test_shadow=sum(test_shadow,2);
            end
            test_shadow=(test_shadow+matrix2_shadow)>0;
            %test_shadow = sum(matrix(:,[ biggest_active diff_order(1:(end-i))]))'>0
            [~,next_shadow] = find_diff_order(test_shadow',matrix(:,diff_order(end))',diff_order(end));
            full_shadow{end+1} = next_shadow{end};
            
        end
    else
        perm_order = [perm_order diff_order(end+1-i)];
        full_shadow{end+1}=full_shadow_test{end-i+1};
        %full_shadow=[full_shadow full_shadow_test(end-i+1)];
        
        new_shadow = (matrix(:,biggest_active)+matrix2_shadow)>0;
        
        [new_diff_order,new_test_shadow]=find_diff_order(new_shadow',matrix(:,diff_order((end-i+2):end))',diff_order((end-i+2):end));
        
        [perm_order,full_shadow]=add_changes(new_diff_order,biggest_active,matrix,perm_order,full_shadow,new_test_shadow,matrix2_shadow);
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
    %full_shadow{(end+1):(end+length(add_shadow))} = add_shadow;
end


end


function actives=find_actives(perm_order_prev)

actives =[];
for i = unique(perm_order_prev)
    actives = [actives mod(sum(perm_order_prev==i)+1,2)];
end

actives=find(actives);
end
