function [flatStrArray] = write(prims, figureBaseStr);
%
% [flatStrArray] = write( prims, figureBaseStr );
%
% Writes the values in quad primitive prims to flatStrArray that is used to upper level
% code, unflattened and saved to the model file.
%
% The size of the quad figure is inferred from the quad primitive array p. 
%

global QUAD_FIGURE;
global FIG_TYPENAMES;

figType	= [ figureBaseStr 'type ' FIG_TYPENAMES{QUAD_FIGURE} ];
flatStrArray = [cellstr(figType) ];

nRows	= size(prims,1);
nCols	= size(prims,2);

for row = 1:nRows
	for col = 1:nCols
        atomBaseStr = [figureBaseStr 'primitive[' num2str(row-1) '][' num2str(col-1) '].'];
		p	= prims(row,col);
		if (row==1 | row==nRows | col==1 | col==nCols)
            atomE = [atomBaseStr 'elongation ' num2str(p.elongation, '%2.17f')];            
            flatStrArray = [flatStrArray ; cellstr(atomE)];
		end
		atomQw = [atomBaseStr 'qw ' num2str(p.q(1), '%2.17f')];
		flatStrArray = [flatStrArray ; cellstr(atomQw)];
		atomQx = [atomBaseStr 'qx ' num2str(p.q(2), '%2.17f')];
		flatStrArray = [flatStrArray ; cellstr(atomQx)];
		atomQy = [atomBaseStr 'qy ' num2str(p.q(3), '%2.17f')];
		flatStrArray = [flatStrArray ; cellstr(atomQy)];
		atomQz = [atomBaseStr 'qz ' num2str(p.q(4), '%2.17f')];
		flatStrArray = [flatStrArray ; cellstr(atomQz)];
		atomR = [atomBaseStr 'r ' num2str(p.r, '%2.17f')];
		flatStrArray = [flatStrArray ; cellstr(atomR)];
		atomT = [atomBaseStr 'theta ' num2str(p.theta*180/pi, '%2.17f')];
		flatStrArray = [flatStrArray ; cellstr(atomT)];

		atomX = [atomBaseStr 'x ' num2str(p.pos(1), '%2.17f')];
		flatStrArray = [flatStrArray ; cellstr(atomX)];
		atomY = [atomBaseStr 'y ' num2str(p.pos(2), '%2.17f')];
		flatStrArray = [flatStrArray ; cellstr(atomY)];
		atomZ = [atomBaseStr 'z ' num2str(p.pos(3), '%2.17f')]; 
		flatStrArray = [flatStrArray ; cellstr(atomZ)];
	end
end

return;
