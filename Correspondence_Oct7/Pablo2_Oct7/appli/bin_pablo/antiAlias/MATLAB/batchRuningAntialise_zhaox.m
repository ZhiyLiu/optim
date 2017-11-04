ddm_dir = 'C:\Pablo2_age\0109\antiAlias\20mhd-ddm\';
movement = [3 3];
files =dir(ddm_dir);
for i=1:size(files,1)
    [suffix outPrefix]= getString_zhaox(files(i).name);
    if strcmp(suffix,'mhd')
        antiAliasWrapper( strcat(ddm_dir, files(i).name), outPrefix, movement );
    end
end