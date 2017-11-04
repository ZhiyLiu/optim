function logQSymPrim = Log(qSymPrim)

if( sum([qSymPrim.inTangentSpace]) ~= 0 )
	qSymPrim
	error( 'symmetric-space Tube Primitive is not in manifold space.' );
end

logQSymPrim = qSymPrim;

for i=1:numel(qSymPrim)
    logQSymPrim(i) = QuadSymPrimitive( qSymPrim.pos, ...
        log(qSymPrim.r), log(qSymPrim.elongation), SphereLog(qSymPrim.Up1),...
        SphereLog(qSymPrim.Um1), true);
end
