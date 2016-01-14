function test_error_save

myfunc
end


function myfunc
for i = 1:100
    k = 3+i;
    if i == 54
        save(strcat('Output/ErrorRuns/error_',num2str(now),'.mat'))
        error('complete')
        
    end
end
end