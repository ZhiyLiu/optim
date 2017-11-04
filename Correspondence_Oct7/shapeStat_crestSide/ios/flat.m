% ==============================================================
% Convert a Pablo registry file to unnested, key-value 
% (flat) text.
% --------------------------------------------------------------
function arrayChar = flat(m3dFilename)
    [pathstr,name,ext] = fileparts(m3dFilename);
    outFilename = 'Flat.txt';
    %outFilename = [ pathstr(1:length(pathstr)-7) '\Flat.txt'];
        
    if (ispc)
    	command = ['flat_d.exe ' m3dFilename ' ' outFilename];
        searchPaths = path;
        while (1)
            [token, searchPaths] = strtok(searchPaths, ';');
            [status, msg] = dos([token '\' command]);
            if(status == 0 || isempty(searchPaths))
                break;
            end
        end
    elseif (isunix)
    	command = ['linflat ' m3dFilename ' ' outFilename];
        searchPaths = path;
        while (1)
            [token, searchPaths] = strtok(searchPaths, ':');
            [status, msg] = system([token '/' command]);
            if(status == 0 || isempty(searchPaths))
                break;
            end
        end
    else
    	disp('Error: This is not a PC or Unix compatible machine!');
        disp('Residue statistics computation not complete!');
        return;
    end
        
    %arrayChar = strvcat(importdata(outFilename, '\n'));
    fid = fopen(outFilename);
    cellArray = textscan(fid, '%s', 'delimiter', '\n');
    arrayChar = char(cellArray{1});
    fclose(fid);
    
    delete 'Flat.txt';
return;
