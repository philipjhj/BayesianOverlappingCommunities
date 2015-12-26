clc
[a,b,c,d]=k_flip2(Q,Z(1,:),4);
a'
b
c{:}
%%

Q'


%% print shadow from none selected elements from Z
nZ = Z(1,:);
nZ(1,a)=0;
sum(Q(:,[find(nZ)]),2)'>0
Z(1,:)
%%
k=4
max_reset=(k-1)*2+1;
reset_points=[1];
for i = 1:k-1%(max_reset-1)/2
    reset_points = [reset_points 1+reset_points(2^(i-1))*2  reset_points];
end
reset_points

%%
clc
[a,b,c,d]=k_flip2(Z,Q(1,:),4)
a

%%

Z'

%% print shadow from none selected elements from Z
nZ = Q(1,:);
nZ(1,a)=0;
sum(Z(:,[find(nZ)]),2)'>0
a'








%%




