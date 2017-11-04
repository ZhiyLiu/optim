function [flatStrArray] = write(prims, figureBaseStr);
%
% [flatStrArray] = write( prims, figureBaseStr );
%
% Writes the values in quad primitive "prims" to flatStrArray that is used
% in upper level code, unflattened and saved to the model file.
%
% The size of the quad figure is inferred from the quad primitive array p. 
%

global QUAD_FIGURE;
global FIG_TYPENAMES;

% dibyendu -  the program cannot find FIG_TYPENAMES and QUAD_FIGURE
% which are declared in the calculatePGA.m file. I am hard-coding the figType here
%
% figType = [ figureBaseStr 'type ' FIG_TYPENAMES{QUAD_FIGURE} ];

figType = [ figureBaseStr 'type ' 'QuadFigure' ];

flatStrArray = [cellstr(figType) ];

nRows   = size(prims,1);
nCols   = size(prims,2);

for row = 1:nRows
    for col = 1:nCols
        atomBaseStr = [figureBaseStr 'primitive[' num2str(row-1) '][' num2str(col-1) '].'];
        p   = prims(row,col);
        
        % dibyendu - this was previously wrong 
        % different strings are reqd for std atoms and end atoms
        % atomType = 0 ==> std atom
        % atomType = 1 ==> end atom

        if( ( row == 1 ) || ( row == nRows ) || ( col == 1 ) || ( col == nCols ) )
            % atomType = 1 ;
            nSpokes = 3 ;
        else
            % atomType = 0 ;
            nSpokes = 2 ;
        end        
        
        for n = 1:nSpokes,
            nm = n - 1;
            atomR = [atomBaseStr 'r[' num2str(nm) '] ' num2str(p.r(n), '%2.17f')];
            flatStrArray = [flatStrArray ; cellstr(atomR)];
            atomUx = [atomBaseStr 'ux[' num2str(nm) '] ' num2str(p.U(1, n), '%2.17f')];
            flatStrArray = [flatStrArray ; cellstr(atomUx)];
            atomUy = [atomBaseStr 'uy[' num2str(nm) '] ' num2str(p.U(2, n), '%2.17f')];
            flatStrArray = [flatStrArray ; cellstr(atomUy)];
            atomUz = [atomBaseStr 'uz[' num2str(nm) '] ' num2str(p.U(3, n), '%2.17f')];
            flatStrArray = [flatStrArray ; cellstr(atomUz)];
        end
        atomX = [atomBaseStr 'x ' num2str(p.pos(1), '%2.17f')];
        flatStrArray = [flatStrArray ; cellstr(atomX)];
        atomY = [atomBaseStr 'y ' num2str(p.pos(2), '%2.17f')];
        flatStrArray = [flatStrArray ; cellstr(atomY)];
        atomZ = [atomBaseStr 'z ' num2str(p.pos(3), '%2.17f')]; 
        flatStrArray = [flatStrArray ; cellstr(atomZ)];
    end
end

return;
