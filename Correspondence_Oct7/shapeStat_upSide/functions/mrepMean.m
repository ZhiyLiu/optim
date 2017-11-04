% Function [M] = mrepMean(Mlist, header)
%
%   mrepMean(Mlist) computes the mean mreps in Mlist.
%
% INPUT : Mlist - cell array of combined(tube, quad) atoms
% OUTPUT : M - cell array of combined(tube, quad) atoms

function M = mrepMean(Mlist, header)

global TUBE_FIGURE;

[nSamps, nAtoms] = size(Mlist);
symRep = 0;
lieRep = 0;
sRep = 0;
if (isa(Mlist{1}, 'QuadSymPrimitive') ||  isa(Mlist{1}, 'TubeSymPrimitive'))
    symRep = 1;
elseif (isa(Mlist{1}, 'QuadPrimitive') ||  isa(Mlist{1}, 'TubePrimitive'))
    lieRep = 1;
elseif isa(Mlist{1}, 'SrepQuadPrimitive') % s-rep tubes not yet supported at UNC
    sRep = 1;
else
    error('Unsupported representation.');
end

for i = 1:nAtoms
    if lieRep
        SM = atomMean( convert2Sym([ Mlist{:, i} ] )); % [ Mlist{:, i} ] row array of primitives of the same type
        M(i) = { convert2Lie(SM) };
    else % symRep or s-rep
        M(i) = { atomMean( [ Mlist{:, i} ] ) };
    end
end

% ===================================================
% if tube figures, align tube figures
% ===================================================
if( lieRep == 1 )
    atomCount	= 1;
    for i = 1:sum(header.nFigs)
        if( header.figTypes(i) == TUBE_FIGURE )
            alignedAtoms = alignTube( ...
                [M{atomCount:atomCount+header.nFigAtoms(i)-1}], ...
                reshape([Mlist{:,atomCount:atomCount+header.nFigAtoms(i)-1}], ...
                [ nSamps, header.nFigAtoms(i)]) );
            for c = 1:size(alignedAtoms,2)
                M(atomCount+c-1) = {alignedAtoms(c)};
            end
        end
        atomCount = atomCount + header.nFigAtoms(i);
    end
end
