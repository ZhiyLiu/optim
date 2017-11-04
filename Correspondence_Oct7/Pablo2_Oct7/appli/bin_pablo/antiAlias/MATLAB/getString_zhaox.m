function [insuffix outPrefix] = getString_zhaox(filename)
filelen=length(filename);
for k=1:filelen
    if filename(k)=='.';
        dotnum=k;
    end
    k=k+1;
end
insuffix= filename(dotnum + 1:end);
outPrefix = filename(1:dotnum-5);
