% ==============================================================
% Function valStr = parsePrimLines(keyStr, modelArray, nFigs, nCols, nRows)
%
%   Input
%       modelArray: a character array of flatted m3d model file
%       nFigs, nCols, nRow: column of integer
%   Output: values in lines corresponding to the key
% 
% NOTE: XXX Need modification for multi-figure objects
% --------------------------------------------------------------

function prims  = parsePrimLines(modelArray, nFigs, nCols, nRows)

% Sample format of prim lines:
%    model.figure[0].primitive[0][0].elongation 1.1875297062069081
%    model.figure[2].primitive[2][6].qx 0.88279279228808361
%

useSrepFormat = (~isempty(strmatch('model.figure[0].primitive[0][0].ux[0]',...
      modelArray)));

numFigs = sum(nFigs);
atomCount = 1;
prims = {};
for figId = 0:(numFigs-1)
	% Check the figure type
	figKeyStr = ['model.figure[' int2str(figId) '].'];    
	figType = findVal([figKeyStr 'type'], modelArray, ' %s');
	
	if (strcmp(figType, 'QuadFigure'))
		for row = 0:(nRows(figId+1)-1)
			rowPrimKeyStr = [figKeyStr 'primitive[' int2str(row) ']'];

			for col = 0:(nCols(figId+1)-1)
				primKeyStr = [rowPrimKeyStr '[' int2str(col) '].'];
                if useSrepFormat
                    quadPrim = SrepQuadPrimitive(primKeyStr, modelArray);
                else
                    quadPrim = QuadPrimitive(primKeyStr, modelArray);
                end

				prims(atomCount) = {quadPrim};

				atomCount = atomCount + 1;
			end
		end
	elseif (strcmp(figType, 'TubeFigure'))
        numberOfSpokes = findVal([figKeyStr 'numberOfSpokes'], modelArray, '%d');
		 for col = 0:(nCols(figId+1)-1)
			primKeyStr = [figKeyStr 'primitive[' int2str(col) '].'];
			tubePrim = TubePrimitive( primKeyStr, modelArray, numberOfSpokes );
			prims(atomCount) = {tubePrim};
			atomCount = atomCount + 1;
		 end
	else
		error('Cannot determine the figure type (QuadFigure or TubeFigure)');
	end
end
return;
