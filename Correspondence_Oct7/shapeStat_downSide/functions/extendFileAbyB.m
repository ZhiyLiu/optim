
function extendFileAbyB(FileA,FileB)

fid_b=fopen(FileB,'r');
fid_a= fopen(FileA, 'at');
b=fread(fid_b);
fwrite(fid_a,b);%write to the filename
fclose(fid_a);
fclose(fid_b);