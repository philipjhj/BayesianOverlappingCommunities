function writeConceptsToFile(cluster_matrix,names,filename,direction,mode)

f = fopen(filename, 'w');

if nargin < 5
    mode = 0;
end

if direction == 'r' %rows
    for i = 1:size(cluster_matrix,2)
        fprintf(f,'%s,',names{cluster_matrix(:,i)==1});
        fprintf(f, '\n');
    end
elseif direction == 'c'
    [row,col]=find(cluster_matrix);
    
    fprintf(f,'Concept %d,',1:size(cluster_matrix,2));
    fprintf(f,'\n');
    
    if mode == 0
        for i = 1:length(unique(row))-1
            for j = 1:size(cluster_matrix,2);
                if sum(col==j)>=i
                    current_row = row(col==j);
                    next_elem = find(current_row,i);
                    next_elem = next_elem(i);
                    fprintf(f,'%s,',names{current_row(next_elem)});
                else
                    fprintf(f,',');
                end
            end
            fprintf(f,'\n');
        end
    else
        for i = 1:size(cluster_matrix,2)
            strings = names(cluster_matrix(:,i)==1);
            japenese=0; norwegian=0; swedish=0;
            for j = 1:length(strings)
         
                string_name = strings{j};
                if string_name(1) == 'j'
                    japenese = japenese+1;
                elseif string_name(1) == 'n'
                    norwegian = norwegian+1;
                elseif string_name(1) == 's'
                    swedish = swedish+1;
                end
            end
            fprintf(f,'%d Japanese \n',japenese);
            fprintf(f,'%d Norwegian \n',norwegian);
            fprintf(f,'%d Swedish \n',swedish);
        end
    end
else
    disp('choose direction with ''r'' (rows) or ''c'' (columns)')
    
end

fclose(f);
end