function [flatStrArray] = write(prims, figureBaseStr);
%
% [flatStrArray] = write( prims, figureBaseStr );
%
% Writes the values in tube primitive array prims to flatStrArray that is used to upper level
% code, unflattened and saved to the model file.
%
% numCols is inferred from the size of the passed array and numRows is always
% 1 for a tube figure.
%

global QUAD_FIGURE;
QUAD_FIGURE = 1;
global TUBE_FIGURE;
TUBE_FIGURE = 2;
global FIG_TYPENAMES;
FIG_TYPENAMES   = { 'QuadFigure', 'TubeFigure' };
figType	= [ figureBaseStr 'type ' FIG_TYPENAMES{TUBE_FIGURE} ];
flatStrArray = [cellstr(figType) ];

nCols	= prod(size(prims));

for col = 1:nCols
	atomBaseStr = [figureBaseStr 'primitive[' num2str(col-1) '].'];
	p	= prims(col);
	if (col==1 | col==nCols)
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

	atomBaseAtom = [atomBaseStr 'baseAtom ' num2str(p.baseAtom, '%1.0f')];
	flatStrArray = [flatStrArray ; cellstr(atomBaseAtom) ];

	atomX = [atomBaseStr 'x ' num2str(p.pos(1), '%2.17f')];
	flatStrArray = [flatStrArray ; cellstr(atomX)];
	atomY = [atomBaseStr 'y ' num2str(p.pos(2), '%2.17f')];
	flatStrArray = [flatStrArray ; cellstr(atomY)];
	atomZ = [atomBaseStr 'z ' num2str(p.pos(3), '%2.17f')]; 
	flatStrArray = [flatStrArray ; cellstr(atomZ)];

	for i=1:numel(p.dr)
		atomDr = [atomBaseStr 'dr[' int2str(i-1) '] ' num2str(p.dr(i), '%2.17f')];
		flatStrArray	= [flatStrArray ; cellstr(atomDr) ];
	end
end

return;
