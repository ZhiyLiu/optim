function [img] = newImage( image, channels, origin, spacing )
%
% function [img] = newImage( image, channels, origin, spacing )
%   Creates a new img structure with the passed in parameters
% and returns it. The structure has the following fields.
%	I					The multi-channel passed in image
%	channels			The number of channels in the image, the last
%						dimension should be of size channels.
%	origin				The passed in origin
%	spacing				The passed in spacing
%	modelToImageScale	The scale for converting from model co-ordinates
%						([0,1]^3) to image co-ordinates
%	modelToImage		A handy function to convert a point from model to
%						image co-ordinates.
%	imageToModel, worldToImage, imageToWorld    other conversion functions.
%	lookupImage			A function that returns an interpolated value from
%						the image at the passed in image co-ordinate.
%	lookup				A faster lookup function ... assumes point to be
%						an integer.
%	fullLookup			A lookup function that takes both an integer point
%						and the channel to lookup into.
% Input:
%   image       A cell array of the actual images.
%	channels	The number of channels in the image
%   origin      The image origin
%   spacing     The voxel size (in model coordinates.)
% Output:
%   img         A structure containing the image and meta information.
%

padding = 0;
dims = size(image);
if( channels > 1 )
	dims = dims(1:end-1);
	if( size(image, ndims(image) ) ~= channels )
		error('newImage: Image does not agree with number of channels');
	end
end
dims = dims - 2*padding;
if (numel(dims) == 3)
	% flip image if any spacing is negative.
	if( spacing(1) < 0 )
		image = image(end:-1:1,:,:);
	end
	if( spacing(2) < 0 )
		image = image(:,end:-1:1,:);
	end
	if( spacing(3) < 0 )
		image = image(:,:,end:-1:1);
	end
	spacing = abs(spacing);
end

extent = dims.*spacing;
modelToImageScale = dims .* max(extent)./extent;

modelToWorld = @(pt) bsxfun(@plus, bsxfun(@times, pt, spacing.*dims), origin - 0.5*spacing);
worldToModel = @(pt) bsxfun(@rdivide, bsxfun(@minus, pt, origin - 0.5*spacing), spacing.*dims);

% the model <-> image coord converstion routines with the 0.5 offset.
% replaced by bsxfun variants for speed.
% modelToImage = @(pt) pt .* repmat(modelToImageScale, [size(pt,1) 1]) - 0.5  + 1 + padding;
% imageToModel = @(pt) (pt - padding - 1 + 0.5)./repmat(modelToImageScale, [size(pt,1) 1]);
% 
% worldToImage = @(pt) (pt-repmat(origin,[size(pt,1) 1]))./repmat(spacing, [size(pt,1) 1]) + 1 + padding;
% imageToWorld = @(pt) (pt - padding - 1).*repmat(spacing, [size(pt,1) 1]) + repmat(origin, [size(pt,1) 1]);

modelToImage = @(pt) bsxfun(@times, pt, modelToImageScale) + (-0.5  + 1 + padding);
imageToModel = @(pt) bsxfun(@rdivide, pt + (-padding - 1 + 0.5), modelToImageScale);

worldToImage = @(pt) bsxfun(@rdivide, bsxfun(@minus, pt, origin), spacing) + (1 + padding);
imageToWorld = @(pt) bsxfun(@plus, bsxfun(@times, pt + (-padding - 1), spacing), origin);

switch(numel(dims))
	case 2
		lookupImage = @lookupImage2;
		fastLookup	= @fastLookup2;
		fullLookup	= @fullLookup2;
	case 3
		lookupImage = @lookupImage3;
		fastLookup	= @fastLookup3;
		fullLookup	= @fullLookup3;
	case 4
		lookupImage = @lookupImage4;
		fastLookup	= @fastLookup4;
		fullLookup	= @fullLookup4;
	otherwise
		warning('image:DimsPartiallySupported', '%d-dimensional images do not have a valid lookupImage function.', numel(dims));
		lookupImage	= [];
		fastLookup	= [];
		fullLookup	= [];
end

img = struct('I', image, 'channels', channels, 'origin', origin, ...
	'spacing', spacing, ...
    'modelToImage', modelToImage, 'imageToModel', imageToModel, ...
    'worldToImage', worldToImage, 'imageToWorld', imageToWorld, ...
    'worldToModel', worldToModel, 'modelToWorld', modelToWorld, ...
	'lookupImage', lookupImage, 'lookup', fastLookup, ...
	'fullLookup', fullLookup, ...
    'modelToImageScale', modelToImageScale);
end

function val = lookupImage2(img, pt)
	%
	% function val = lookupImage2(img, pt)
	%	An internal function that returns an interpolated value from a 2D
	% image at the given point in image co-ordinates. The interpolation is
	% linear and trying to extrapolate results in NaNs.
	%
	% Input:
	%	img	A 2D array representing a single channel through the
	%		multi-channel 2D image.
	%	pt	The image co-ordinates of the point to be looked up.
	%		pt should be 1x2 or 2x1.
	% Output:
	%	val	The interpolated value of the image at the point.
	%
	ipt = floor(pt);
	pt = pt - ipt + 1;
	val = interpn( ...
		img( ipt(1):ipt(1)+1, ipt(2):ipt(2)+1 ), ...
		pt(1), pt(2) );
end

function val = lookupImage3(img, pt)
	%
	% function val = lookupImage3(img, pt)
	%	An internal function that returns an interpolated value from a 3D
	% image at the given point in image co-ordinates. The interpolation is
	% linear and trying to extrapolate results in NaNs.
	%
	% Input:
	%	img	A 3D array representing a single channel through the
	%		multi-channel 3D image.
	%	pt	The image co-ordinates of the point to be looked up.
	%		pt should be 1x3 or 3x1.
	% Output:
	%	val	The interpolated value of the image at the point.
	%
	ipt = floor(pt);
	pt = pt - ipt + 1;
	val = interpn( ...
		img( ipt(1):ipt(1)+1, ipt(2):ipt(2)+1, ipt(3):ipt(3)+1 ), ...
		pt(1), pt(2), pt(3) );
end

function val = lookupImage4(img, pt)
	%
	% function val = lookupImage4(img, pt)
	%	An internal function that returns an interpolated value from a 4D
	% image at the given point in image co-ordinates. The interpolation is
	% linear and trying to extrapolate results in NaNs.
	%
	% Input:
	%	img	A 4D array representing a single channel through the
	%		multi-channel 4D image.
	%	pt	The image co-ordinates of the point to be looked up.
	%		pt should be 1x4 or 4x1.
	% Output:
	%	val	The interpolated value of the image at the point.
	%
	ipt = floor(pt);
	pt = pt - ipt + 1;
	% some special cases
	if( abs(pt(4) - 1) <= 1e-4 )
		pt = pt + 1;
		if( ipt(4) == 1 || ipt(4) == size(img,4) )
			val = interpn( ...
				img( ipt(1)-1:ipt(1)+1, ipt(2)-1:ipt(2)+1, ipt(3)-1:ipt(3)+1, ipt(4) ), ...
				pt(1), pt(2), pt(3), 'cubic' );
		else
			val = interpn( ...
				img( ipt(1)-1:ipt(1)+1, ipt(2)-1:ipt(2)+1, ipt(3)-1:ipt(3)+1, ipt(4)-1:ipt(4)+1 ), ...
				pt(1), pt(2), pt(3), pt(4), 'cubic' );
		end
	else
		val = interpn( ...
			img( ipt(1):ipt(1)+1, ipt(2):ipt(2)+1, ipt(3):ipt(3)+1, ipt(4):ipt(4)+1 ), ...
			pt(1), pt(2), pt(3), pt(4) );
	end
end


%
% These series of functions are similar to the above functions, except that
% there is no interpolation done and pt is assumed to contain integer
% values.
%
function val = fastLookup2(img, pt)
	val = squeeze(img( pt(1), pt(2), : ));
end

function val = fastLookup3(img, pt)
	val = squeeze(img( pt(1), pt(2), pt(3), : ));
end

function val = fastLookup4(img, pt)
	val = squeeze(img( pt(1), pt(2), pt(3), pt(4), : ));
end

function val = fullLookup2(img, pt, channel)
	val = img( pt(1), pt(2), channel );
end

function val = fullLookup3(img, pt, channel)
	val = img( pt(1), pt(2), pt(3), channel );
end

function val = fullLookup4(img, pt, channel)
	val = img( pt(1), pt(2), pt(3), pt(4), channel );
end
