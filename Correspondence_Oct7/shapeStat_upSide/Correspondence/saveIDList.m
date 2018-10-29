function saveIDList(filename, idList)
    fid = fopen(filename, 'w+');
    if(fid < 0)
        return;
    end
    
    for i = 1:size(idList)
        fprintf(fid, '%d\n', idList(i));
    end
    
    fclose(fid);
end