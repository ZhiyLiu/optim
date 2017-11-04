% Function M = makeDataMtx(atomMatrix, varargin)
%
% 1) WHEN varargin is empty,
%
%   makeDataMtx(atomMatrix) converts cell array of m-rep primitives into a
%   double matrix of atoms.
%
% INPUT : atomMatrix - cell array of combined(tube, quad) atoms
% OUTPUT : M - double matrix of combined(tube, quad) atoms
%
% 2) WHEN NOT empty,
%
%   makeDataMtx(atomMatrix) converts a double matrix of atoms into cell
%   array of m-rep primitives.
%
% INPUT : atomMatrix - double array of combined(tube, quad) atoms in
%                      TANGENT SPACE   [nSamps X numVars(dims*nAtoms)]
%            
%         varargin{1} - int array of indices whether each atom is in
%                       tube(0) or quad(1)  [1 X nAtoms]
% OUTPUT : M - cell array of combined(tube, quad) atoms
%
% NOTE!!!!!!  DON'T FORGET atoms are in TANGENT SPACE!!!
%
%
% Rohit: For tubes, should we include the spoke variation in here, if
% we do then the figure level PGA will also move the spokes around.
% If not, then no sweat.
%

function [M, varargout] = makeDataMtx(atomMatrix, varargin)

if (nargin == 1)

    [row, col] = size(atomMatrix);

    M = [];

    for r = 1:row
        aRow = [];
        for c = 1:col
            aRow = [aRow, cVec([ atomMatrix{r,c} ])];   
            % col X ((tubeSym)9 or (quadSym)11)
        end
        [nAtoms, aDims] = size(aRow);        % col == nAtoms
        M = [M ; reshape(aRow', [1, nAtoms*aDims])] ;
    end

    if (nargout > 1 )

        atomIndicator = zeros(1,col);

        for c = 1:col
            if (isa(atomMatrix{1,c}, 'QuadSymPrimitive'))
                atomIndicator(c) = 1;
            elseif (isa(atomMatrix{1,c}, 'TubeSymPrimitive'))
                atomIndicator(c) = 0;
            else
                disp('   Error: Not valid m-rep primitive.');
                disp(' ');
                return;
            end
        end
        
        varargout{1} = atomIndicator;
    end

else
    
    [row, col] = size(atomMatrix);
    atomIndicator = varargin{1};
    nAtoms = length(varargin{1});

    
    for r = 1:row
        
        icol = 1;
        
        for c = 1:nAtoms
            
            if ( atomIndicator(c) == 1 )  % Quad (9 parameters)
                M{r, c} = QuadSymPrimitive( atomMatrix(r, icol:icol+2)', atomMatrix(r, icol+3), ...
                                            atomMatrix(r, icol+8), atomMatrix(r, icol+4:icol+5)', ...
                                            atomMatrix(r, icol+6:icol+7)', true);
                icol = icol + 9;
            elseif ( atomIndicator(c) == 0 ) % Tube (8 parameters)
                M{r, c} = TubeSymPrimitive( atomMatrix(r, icol:icol+2)', atomMatrix(r, icol+3), ...
                                            atomMatrix(r,  icol+7), atomMatrix(r, icol+4:icol+5)', ...
                                            atomMatrix(r, icol+6), false, ...
											repmat(0.0, [8,1]), ...	% TODO: Do something about this!!!
											true);
                icol = icol + 8;                                        
            else
                disp('   Error: Not valid value in the 2nd argument.');
                disp(' ');
                return;
            end
                
        end
    end
    
end

return;
