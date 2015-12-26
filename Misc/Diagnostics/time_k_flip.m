load('time_k_flip_200_200_25.mat')

for k = 2:20%size(Z,2)
    for i = 1:size(Z,1)
        tic;
        [~,~,~,~] = k_flip2(Q,Z(i,:),k);
        t_k(k-1,i)=toc;
    end
end


save('time_of_k_flip_200_200_25','t_k')