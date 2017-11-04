function qSymPrim = adjustUnits(qSymPrim, mr, bCommensurate, varargin)

% bCommensurate -> 1 : multiplication, bCommensurate -> 0 : division 
% mr: mean radius 

weightMode = varargin{1};
if (weightMode) 
    w = mr/2;
else
    w = mr;
end
if( qSymPrim.inTangentSpace == 1 )
	if (bCommensurate)
		qSymPrim = QuadSymPrimitive(qSymPrim.pos,...
			mr*qSymPrim.r, w*qSymPrim.elongation, w*qSymPrim.Up1, w*qSymPrim.Um1, true);
	else
		qSymPrim = QuadSymPrimitive(qSymPrim.pos,...
			qSymPrim.r/mr, qSymPrim.elongation/w, qSymPrim.Up1/w, qSymPrim.Um1/w, true);
	end
else
	qSymPrim
	error('QuadSymPrimitive prim should be in tangent space.');
end

return;









