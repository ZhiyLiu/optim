function quadSymPrim = QuadSymPrimitive(varargin)

% quadPrim = QuadPrimitive(pos: 3x1 vector, r, elongation, Up1: 3x1(or 2x1) vector,
%                          Um1: 3x1(or 2x1) vector, inTangentSpace: boolean)

switch nargin 
    case 0
        quadSymPrim = struct('pos', [0.0; 0.0; 0.0], 'r', 1.0, 'elongation', 1.0,...
            'Up1', [1.0; 0.0; 0.0], 'Um1', [1.0; 0.0; 0.0], 'inTangentSpace', false);
        quadSymPrim = class(quadSymPrim, 'QuadSymPrimitive');
    case 1
        if (isa(varargin{1}, 'QuadSymPrimitive'))
            quadSymPrim = varargin{1};
        else 
            error('Wrong argument type');
        end
    case 5
        if (length(varargin{4}) == 2 && length(varargin{4}) == 2)
            quadSymPrim = struct('pos', varargin{1} , 'r', varargin{2}, 'elongation', varargin{3},...
                'Up1', varargin{4}, 'Um1', varargin{5}, 'inTangentSpace', true);
        elseif (length(varargin{4}) == 3 && length(varargin{4}) == 3)
            quadSymPrim = struct('pos', varargin{1} , 'r', varargin{2}, 'elongation', varargin{3},...
                'Up1', varargin{4}, 'Um1', varargin{5}, 'inTangentSpace', false);
        else
            error('Wrong argument');
        end
        quadSymPrim = class(quadSymPrim, 'QuadSymPrimitive');           
    case 6
        quadSymPrim = struct('pos', varargin{1} , 'r', varargin{2}, 'elongation', varargin{3},...
            'Up1', varargin{4}, 'Um1', varargin{5}, 'inTangentSpace', varargin{6});
        quadSymPrim = class(quadSymPrim, 'QuadSymPrimitive');
    otherwise
        error('Wrong number of input arguments')
end
    
