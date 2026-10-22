function json(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>0);clean=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
end
