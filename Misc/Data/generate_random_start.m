A = wvs_japan_norway_sweden_males;
%%

d = 25;

for i = 1:5
Z{i} = randi([0 1],size(A,1),d);
Q{i} = randi([0 1],size(A,2),d);
end

%%

save('wvs_random_intial_Z_Q','Z','Q')
