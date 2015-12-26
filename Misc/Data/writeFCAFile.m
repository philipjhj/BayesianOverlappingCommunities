function writeFCAFile(M,obj,fea,filename)

f = fopen(filename, 'w');

% Write first line with features
fprintf(f,';');
fprintf(f,'%s;',fea{:});
fprintf(f,'\n');

% Write objects and corresponding relations
for i = 1:length(obj)
    fprintf(f, '%s;', obj{i});
    fprintf(f, '%d;', M(i,:));
    fprintf(f, '\n');
end
fclose(f);