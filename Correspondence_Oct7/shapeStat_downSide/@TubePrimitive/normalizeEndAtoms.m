function [prims] = normalizeEndAtoms( prims )
%
%  [normalizedAtoms] = normalizeEndAtoms(prims)
% 
% Takes a 2D array of Tube primitives and normalizes the end-atoms
% The 2D array should belong to the same population and consist of just one tube figure.
%

%
% First, convert all the m-reps to symmetric space.
%
for sample=1:size(prims,1)
	sym(sample,:)	= convert2Sym(prims(sample,:));
end

% atom 0 is end with smaller index
% atom 1 is end with larger index
for atom = 0:1
	if( atom == 0 )
		iult	= 1;
		ipenult	= 2;
		i2penult= 3;
		istep	= +1;
	else
		iult	= size( prims, 2);
		ipenult	= iult-1;
		i2penult= ipenult-1;
		istep	= -1;
	end

	%%
	%% new method
	%%

	%
	% Compute mean arc length.
	%
	temp			= get(sym(:,iult),'pos') - get(sym(:,ipenult),'pos');
	arcLengthDiff	= sqrt(sum(temp.*temp,1));
	meanArcLength	= mean(arcLengthDiff);
	%
	% Now find interpolation/extrapolation parameter for each m-rep.
	%
	t	= repmat(meanArcLength,size(arcLengthDiff)) ./ arcLengthDiff;
	%
	% Find the difference of the ultimate and penultimate atoms in
	% symmetric space (diff is in coords of penult atom)
	%
	penultUltDiff	= Log( sym(:,iult) - sym(:,ipenult) );
	%
	% Now, replace the end-atom with the penultimate + interpolated
	% version of end-atom.
	%
	sym(:,iult)		= Exp( Log(sym(:,ipenult)) + penultUltDiff .* t );

	%%
	%% OLD method
	%%
% 		%
% 		% Extract the two end-atoms from one end
% 		%
% 		endhead	= header;
% 		endraw	= rawDataStruct;
% 		endhead.nAtoms	= 2;
% 		endraw.atoms	= endraw.atoms(:,:,ipenult:istep:i2penult);
% 
% 		%
% 		% Align the penultimate atoms in trans and rotation.
% 		% alignMode 2 is trans+rot.
% 		alignedDataStruct = alignModels(2, endhead, endraw);
% 		%
% 		% Now align the ultimate atoms (aligning the whole m-rep in the
% 		% process) with the alignment computed in the previous step.
% 		%
% 		ultimate	= alignedDataStruct;
% 		ultimate.atoms	= rawDataStruct.atoms;
% 		for samp = 1:header.nSamps
% 			ultimate.atoms(:,samp,:)	= reshape( ...
% 				applyMrepTransf( ultimate.alignTransf(:,samp), squeeze(ultimate.atoms(:,samp,:)) ), ...
% 				[size(rawDataStruct.atoms,1), 1, header.nAtoms] );
% 		end
% 
% 		% now convert all the m-reps to symmetric space.
% 		ultimate.atoms = convert2Sym(ultimate.atoms, header.figTypes(1) );
% 
% 		%
% 		% Find the mean of the difference of the ultimate and penultimate atoms
% 		%
% 		penultUltDiff	= atomDiff( ultimate.atoms(:,:,ipenult), ultimate.atoms(:,:,iult), header.figTypes(1) );
% 		MA	= atomMeanSym( penultUltDiff, header.figTypes(1) );
% 
% 		%
% 		% Now replace the end-atom with the penultimate atom + mean
% 		% diff of (ultimate atom, penultimate atom) (diff is in coords
% 		% of penultimate atom)
% 		%
% 		ultimate.atoms(:,:,iult)	= atomMult( ultimate.atoms(:,:,ipenult), ...
% 			reshape( repmat(MA, header.nSamps, 1), size( ultimate.atoms(:,:,ipenult)) ), ...
% 			header.figTypes(1) );
% 		%
% 		% Transform back to lie group space and save in header.
% 		%
% 		rawDataStruct.atoms	= convert2Atom( ultimate.atoms, header.figTypes(1) );
end

%
% Transform back to lie group space
%
for sample=1:size(prims,1)
	tube	= convert2Lie(sym(sample,:));
	prims(sample,:)	= alignTube(tube, prims(sample,:));
end


