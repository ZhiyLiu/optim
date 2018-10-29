% read linking structure and parse into a column vector( 68*3-dimensional, then return this
% vector
% input: single s-rep file path
function f = loadLinkingStructure(file_path, file_name)
    data_file = fullfile(file_path, file_name);
    
    fid = fopen(data_file, 'r');
    f = [];
    nObjs = 0;

    objId = -1;
    nLinks = 0;
    
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        %disp(tline);
%         [s, e] = regexp(tline, 'numOfObj = ');
%       
%         nObjs = tline(e+1);
        res = regexp(tline, 'numOfObj = (?<n>\d+);', 'names');
        if size(res) ~= 0
           nObjs = res.n;
        end
        
        res = regexp(tline, 'figure[(?<id>\d+)] {', 'names');
        if size(res) ~= 0
            objId = res.id;
        end

        % the number of links of this object id
        res = regexp(tline, 'numOfLinks = (?<nLinks>\d+);', 'names');
        if size(res) ~= 0
            nLinks = res.nLinks;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% position value of this link  %%%%%%%%%%%%
        res = regexp(tline, '\<x = (?<val>-?\d+\.?(\d+)?(e-)?(\d+)?);', 'names');
        if size(res) ~= 0
            f = [f; str2num(res.val)];
        end
        
        res = regexp(tline, '\<y = (?<val>-?\d+\.?(\d+)?(e-)?(\d+)?);', 'names');
        if size(res) ~= 0
            f = [f; str2num(res.val)];
        end
        
        res = regexp(tline, '\<z = (?<val>-?\d+\.?(\d+)?(e-)?(\d+)?);', 'names');
        if size(res) ~= 0
            f = [f; str2num(res.val)];
        end
        
        %%%%%%%%%%%%%%%%% spoke direction %%%%%%%%%%%%%%%%%%%%%%%%%%
        res = regexp(tline, '\<ux = (?<val>-?\d+\.?(\d+)?(e-)?(\d+)?);', 'names');
        if size(res) ~= 0
            f = [f; str2num(res.val)];
        end
        
        res = regexp(tline, '\<uy = (?<val>-?\d+\.?(\d+)?(e-)?(\d+)?);', 'names');
        if size(res) ~= 0
            f = [f; str2num(res.val)];
        end
        
        res = regexp(tline, '\<uz = (?<val>-?\d+\.?(\d+)?(e-)?(\d+)?);', 'names');
        if size(res) ~= 0
            f = [f; str2num(res.val)];
        end
        
        %%%%%%%%%%%%%% spoke length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        res = regexp(tline, '\<r = (?<val>-?\d+\.?(\d+)?(e-)?(\d+)?);', 'names');
        if size(res) ~= 0
            f = [f; str2num(res.val)];
        end
        
        %%%%%%%%%%%%%%% link length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        res = regexp(tline, '\<linkLength = (?<val>-?\d+\.?(\d+)?(e-)?(\d+)?);', 'names');
        if size(res) ~= 0
            f = [f; str2num(res.val)];
        end
        
%         %%%%%%%%%%%%%%% Link point z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         res = regexp(tline, '\<linkPoint_x = (?<val>-?\d+\.?(\d+)?);', 'names');
%         if size(res) ~= 0
%             f = [f; str2num(res.val)];
%         end
%         
%         res = regexp(tline, '\<linkPoint_y = (?<val>-?\d+\.?(\d+)?);', 'names');
%         if size(res) ~= 0
%             f = [f; str2num(res.val)];
%         end
%         
%         res = regexp(tline, '\<linkPoint_z = (?<val>-?\d+\.?(\d+)?);', 'names');
%         if size(res) ~= 0
%             f = [f; str2num(res.val)];
%         end
        %%%%%%%%%%%%%%% Link vector from z' to z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        res = regexp(tline, '\<linkVector_x = (?<val>-?\d+\.?(\d+)?(e-)?(\d+)?);', 'names');
        if size(res) ~= 0
            f = [f; str2num(res.val)];
        end
        
        res = regexp(tline, '\<linkVector_y = (?<val>-?\d+\.?(\d+)?(e-)?(\d+)?);', 'names');
        if size(res) ~= 0
            f = [f; str2num(res.val)];
        end
        
        res = regexp(tline, '\<linkVector_z = (?<val>-?\d+\.?(\d+)?(e-)?(\d+)?);', 'names');
        if size(res) ~= 0
            f = [f;0];%[f; str2num(res.val)];
        end
        
        %%%%%%%%%%%%%%%%%%%% Additional link information %%%%%%%%%%%%%
        res = regexp(tline, '\<linkTo = (?<val>-?\d+\.?(\d+)?);', 'names');
        if size(res) ~= 0
            f = [f; str2num(res.val)];
        end
        
        % parse each figure (obj) in the configuration
%         for j = 0: nObjs-1
%             start_figure = strcat('figure[', j, ']');
%             if size(strfind(tline, start_figure)) > 0
%                 % start to parse this object
%                 
%             end
%         end

    end
    

    fclose(fid);
    
end