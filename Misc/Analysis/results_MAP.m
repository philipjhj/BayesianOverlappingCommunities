


%files = dir('Results_cluster/Results_old_local/davis_d6*');
%files = dir('Results_cluster/Results_final/davis_d*');
files = dir('Results_cluster/Results_final/davis_d15_*');
%files = dir('Results_cluster/Results_final/AWA*');


MAP = -Inf;
figure
hold on
i=1;
cc = hsv(numel({files.name}));
for file = files'
    load(file.name);
    max_prob = Find_max(logPs);
%    disp(max_prob)
    
    if max_prob>MAP
        MAP = max_prob;
        best_run = file.name;
    end
    
    %plot LogPs
    plot(logPs(:,3),'color',cc(i,:))
    i=i+1;
    
end
disp(best_run)
load(best_run);
[max_prob, sample_max_Z, sample_max_Q] = Find_max(logPs);
disp(max_prob)

Z_best=dZ{sample_max_Z};
Q_best=dQ{sample_max_Q};

figure
[Z_idx,Q_idx,Z_idx_cluster,Q_idx_cluster]=plotZQ(Z_best,Q_best,0,1,women,events);
figure


Z_adjusted = fliplr(Z_best(Z_idx,Z_idx_cluster).*repmat((Z_idx_cluster*1/2+1),size(Z_best,1),1));
Q_adjusted = flipud((Q_best(Q_idx,Q_idx_cluster).*repmat(Q_idx_cluster*1/2+1,size(Q_best,1),1))');

img=[zeros(size(Z_best,2)) Q_adjusted; Z_adjusted A(Z_idx,Q_idx)];

imagesc(img) 
 