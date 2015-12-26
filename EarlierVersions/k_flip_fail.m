function permutations=k_flip(Q,k)
% permutations, permutation_full_shadow

%k=size(Q,2);%4;

[perm_order,~,full_shadow,rand_idx]=find_perm_order(Q,k)


permutations(1,:) = ones(1,k);
for i = 1:length(perm_order)
    permutations(i+1,:) = permutations(i,:);
    permutations(i+1,find(rand_idx==perm_order(i)))=abs(permutations(i+1,find(rand_idx==perm_order(i)))-1);
end

permutations;
end

function [perm_order,current_idx,full_shadow,rand_idx] = find_perm_order(Q,k,perm_order,current_idx,rand_idx,full_shadow)
    
    
    if ~exist('rand_idx','var')
        d=size(Q,2);
        rand_idx = sort(randsample(d,k));
    end
    rand_idx=sort(rand_idx);
    most_active = sum(Q(:,rand_idx),1);
    [~, order] = sort(most_active,'descend');
    flip_order = rand_idx(order);
    
if ~exist('perm_order','var')
    perm_order = zeros(1,2^k-1);
    full_shadow = zeros(1,2^k-1);
    current_idx=1;
end

n_flip = length(flip_order);

used_idx = zeros(1,n_flip);
used_idx(1) = flip_order(1);
test_idx = flip_order(2:end);
n_test = length(test_idx);
diff_order(1) = flip_order(1);

for i = 1:n_test
    [diff_order(i+1),permutation_full_shadow(i)]=most_different((sum(Q(:,used_idx(1:i)),2)>0)', ...
        Q(:,test_idx)',test_idx);
    test_idx(find(test_idx==diff_order(i+1),1)) = [];
    used_idx(i+1) = diff_order(i+1);
end

if n_flip>2
for i = 1:length(diff_order)
    if i < 3
    perm_order(current_idx:(current_idx+i-1))=diff_order((end-i+1):end);
    full_shadow(current_idx:(current_idx+i-1))=permutation_full_shadow((end-i+1):end);
    current_idx=current_idx+i;
    else
        perm_order(current_idx) = diff_order(end+1-i);
        %full_shadow(current_idx)=permutation_full_shadow(end+1-i);
        current_idx=current_idx+1;
%         diff_order((end-i+1):end)
        [perm_order,current_idx,full_shadow] = find_perm_order(Q,k,perm_order,...
            current_idx,diff_order((end-i+2):end),full_shadow);
    end
end
else
    for i = 1:length(flip_order)
    perm_order(current_idx:(current_idx+i-1))=flip_order((end-i+1):end);
    %full_shadow(current_idx:(current_idx+i-1))=full_shadow((current_idx-k):(current_idx-k+i-1));
    current_idx=current_idx+i;
    end
    
end


if length(flip_order)>2
perm_order(current_idx)=flip_order(1);
current_idx=current_idx+1;
[perm_order,~,full_shadow]=find_perm_order(Q,k,perm_order,current_idx,flip_order(2:end),full_shadow);
end


   %perm_order
% if n_flip>2
%     local_idx = current_idx;
% for i = 1:(n_flip-1)
%     perm_order((local_idx-i):(local_idx-1))
%     [perm_order,current_idx]=find_perm_order(Q,k,perm_order,current_idx,...
%     diff_order((local_idx-i):(local_idx-1)));
%     perm_order(current_idx)=perm_order(local_idx-1-i);
%     current_idx=current_idx+1;
% end
%end
end

function [number,full_shadow] = most_different(current,test,test_idx)

reps=size(test,1);
most_shadow = sum((repmat(current,reps,1)-test)==-1,2);

[values,order]=sort(most_shadow,'descend');
number=test_idx(order(1));
if values(1) == 0;
    full_shadow = 1;
else
    full_shadow = 0;
end

end

    