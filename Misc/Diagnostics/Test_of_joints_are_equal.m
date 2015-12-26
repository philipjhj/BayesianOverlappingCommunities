clear all;

load('Results_testing\Gibbs_sampler_testing\Small_data\generated_data01')

suff = suff_stats_class;
p_n(1)=logjoint(A,Z,Q,pms);
[p_opt(1) NZQ1 suff_stat] = logjoint_opt_2(A,Z,Q,pms,suff);

dQ = Q;
dQ(1,1) = abs(dQ(1,1)-1);

p_n(2)=logjoint(A,Z,dQ,pms)
[p_opt(2) NZQ2 suff_stat2] = logjoint_opt_2(A,Z,dQ,pms,suff_stat,'q',[1 1],NZQ1,Q)

dQ2 = Q;
dQ2(1,2) = abs(dQ(1,2)-1);
%clear logjoint_opt
p_n(3)=logjoint(A,Z,dQ2,pms)
[p_opt(3) NZQ3 suff_stat3] = logjoint_opt_2(A,Z,dQ2,pms,suff_stat,'q',[1 2],NZQ1,Q)

% dQ2 = Q;
% dQ2(1,2) = abs(dQ(1,2));
% %clear logjoint_opt
% p_n(3)=logjoint(A,Z,dQ2,pms)
% [p_opt(3) NZQ3] = logjoint_opt(A,Z,dQ2,pms,'q',[1 2],NZQ1,Z,Q)
% %clear logjoint_opt
dZ = Z;
dZ(1,1) = abs(dZ(1,1)-1);
p_n(4)=logjoint(A,dZ,dQ,pms)
[p_opt(4) NZQ4 suff_stat4]=logjoint_opt_2(A,dZ,dQ,pms,suff_stat2,'z',[1 1],NZQ2,Z)

dQ4 = dQ;
dQ4(1,4) = abs(dQ4(1,4)-1);
p_n(5)=logjoint(A,dZ,dQ4,pms)
p_opt(5)=logjoint_opt_2(A,dZ,dQ4,pms,suff_stat4,'q',[1 4],NZQ4,Q)
