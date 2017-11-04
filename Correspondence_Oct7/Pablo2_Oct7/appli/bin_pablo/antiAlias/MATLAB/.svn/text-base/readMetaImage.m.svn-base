function [img] = readMetaImage(filename, headerOnly, verbose)
%
% function [img] = readMetaImage(filename, headerOnly, verbose)
%	Loads 2D and 3D META images and vector fields
% Input
%	filename	The file name of the meta image.
%	headerOnly  An optional argument (default false) to read in only the
%				header.
%	verbose		An optional argument (default false) to display more
%				information.
% Output
%	img			The image in the form of an image structure. See newImage
%				for help on the structure.
%
if nargin < 2
	headerOnly = false;
end
if nargin < 3
	verbose = false;
end

dataFilename='';
dataType = '';
imSize = [0 0 0];
numChannels = 1;
origin = [0 0 0];
spacing = [1 1 1];
headerSize= -1;
swapBytes = false;
% parse header file
fprintf('Loading image: %s...\n',filename);
[key,val] = textread(filename,'%s=%[^\n]');
for i=1:size(key,1)
	switch key{i}
		case 'ObjectType'
			if verbose
				fprintf('  Object Type: %s\n', val{i});
			end
		case 'NDims'
			if verbose
				fprintf('  Number of Dimensions: %s\n', val{i});
			end
		case 'DimSize'
			if verbose
				fprintf('  Size: %s\n', val{i});
			end
			imSize = str2num(val{i});
		case 'ElementType'
			if verbose
				fprintf('  Element Type: %s\n', val{i});
			end
			dataType = decideMETADataType(val{i});
		case 'ElementDataFile'
			if verbose
				fprintf('  DataFile: %s\n', val{i});
			end
			dataFilename = val{i};
		case 'ElementNumberOfChannels'
			if verbose
				fprintf('  Number of Channels: %s\n', val{i});
			end
			numChannels = str2double(val{i});
		case 'Offset'
			if verbose
				fprintf('  Offset: %s\n', val{i});
			end
			origin = str2num(val{i});
		case 'ElementSpacing'
			if verbose
				fprintf('  Spacing: %s\n', val{i});
			end
			spacing = str2num(val{i});
		case 'HeaderSize'
			if verbose
				fprintf('  HeaderSize: %s\n', val{i});
			end
			headerSize = str2double(val{i});
		case 'ElementByteOrderMSB'
			if verbose
				fprintf('  ElementByteOrderMSB: %s\n', val{i});
			end
			swapBytes = strcmpi(val{i}, 'True');
		case 'BinaryDataByteOrderMSB'
			if verbose
				fprintf('  BinaryDataByteOrderMSB: %s\n', val{i});
			end
			swapBytes = strcmpi(val{i}, 'True');
		otherwise
			if verbose
				fprintf('  Unknown Key/Val: %s/%s\n',key{i},val{i});
			end
	end
end

if (headerOnly)
	I = zeros([imSize numChannels],'uint8');
else
	% load image data
	path = fileparts(filename);

	fid = fopen([path filesep dataFilename],'r');
	if (fid==-1)
		error('Can''t read file: %s\n', [path filesep dataFilename]);
	end
	if (headerSize ~= -1)
		fseek(fid, headerSize, 'bof');
	else
		% header size is not specified, read from the end of the file
		% when raw3--> meta, the headersize is not specified, (should be 88 for lung images)
		switch dataType
			case {'float32','float'}
				bytes = 4;
			case {'uint8'}
				bytes=1;
			case{'uint16','int16'}
				bytes =2;
			case{'float64'}
				bytes = 8;
		end

		bytesExpected = prod(imSize)*numChannels*bytes;
		fseek( fid, -bytesExpected, 'eof' );
	end

	% headerSize = ftell(fid);

	[I, count] = fread(fid,prod(imSize)*numChannels,['*' dataType]);
	fclose(fid);
	if( swapBytes)
		I = swapbytes(I);
	end
	I = reshape(I,[numChannels imSize]);
	% This type of image structure needed image data to be aligned such that
	% channels is the last dimension.
	I = shiftdim(I,1);
	% This type of image structure needs image data such that each channel is a
	% separate cell in a cell array.
	% I = num2cell(I, 2:ndims(I) );

	I = double(I);
end

img = newImage(I, numChannels, origin, spacing);
