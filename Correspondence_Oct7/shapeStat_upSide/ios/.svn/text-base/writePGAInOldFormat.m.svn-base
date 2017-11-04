% Writing statistics in Tom's old format

function writePGAInOldFormat(PabloPGAFile, statsFlag, header, globalStruct)

fid = fopen(PabloPGAFile, 'wt');
if (fid == 0)
    error(['Cannot open file ' PabloPGAFile]);
end

pga = globalStruct.stats.PCs(:, 1:globalStruct.projDims);
PGLength = size(pga, 1);

if (statsFlag == 1)
    NUMPARAM = 9;
elseif (statsFlag == 0)
    NUMPARAM = 8;
end

if (header.nAtoms*NUMPARAM ~= PGLength)
    error('Dimension is not matched.');
end

fprintf(fid, '%d\n', globalStruct.projDims);
fprintf(fid, '%d\n', PGLength);


pga = pga(:);
for i = 1:length(pga)
    fprintf(fid, '%12.16f\n', pga(i));
end


fclose(fid);