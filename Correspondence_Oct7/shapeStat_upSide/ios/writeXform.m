% Xform: 8 by number of Xform matrix
%         In each column, the first 4 elements are quarternion, 5th is
%         scale, and the last 3 are translation.
% KeepInv: flag to indicate whether to record Xform (0) or its inverse (1)
% filnames: number of filenames should be the same as the number of columns
%           of Xform matrix

function writeXform(outputDir, filenames, Xform , keepInv)

if (numel(filenames) ~= size(Xform, 2))
    disp('Error in writeXform')
    disp(' ')
    return;
end

for i=1:size(Xform,2)
    Q = Xform(1:4, i);
    r = Xform(5, i);
    x = Xform(6:8, i);
    fid = fopen(fullfile(outputDir, filenames{i}), 'w+');
    if (fid == -1)
        disp(['Error in opening the file: ' fullfile(outputDir, filenames{i})]);
        disp('');
        return;
    end
    
    fprintf(fid, '%s\n', 'model {');
    fprintf(fid, '%s\n', '   figureCount = 0;');
    fprintf(fid, '%s\n', '   name = default;');
    fprintf(fid, '%s\n%s\n%s\n%s\n', '   figureTrees {', '      count = 0;', '   }', '   transformation {');
    
    if (keepInv)
        fprintf(fid, '%s = %1.16f;\n', '      scale', r^(-1));
        invQ = QuatInv(Q);
        fprintf(fid, '%s\n', '      rotation {');
        fprintf(fid, '%s = %1.16f;\n', '         w', invQ(1));
        fprintf(fid, '%s = %1.16f;\n', '         x', invQ(2));
        fprintf(fid, '%s = %1.16f;\n', '         y', invQ(3));
        fprintf(fid, '%s = %1.16f;\n', '         z', invQ(4));
        fprintf(fid, '%s\n', '      }');
        fprintf(fid, '%s\n', '      translation {');
        fprintf(fid, '%s = %1.16f;\n', '         x', -x(1));
        fprintf(fid, '%s = %1.16f;\n', '         y', -x(2));
        fprintf(fid, '%s = %1.16f;\n', '         z', -x(3));
        fprintf(fid, '%s\n', '      }');
    else
        fprintf(fid, '%s = %1.16f;\n', '      scale', r);
        fprintf(fid, '%s\n', '      rotation {');
        fprintf(fid, '%s = %1.16f;\n', '         w', Q(1));
        fprintf(fid, '%s = %1.16f;\n', '         x', Q(2));
        fprintf(fid, '%s = %1.16f;\n', '         y', Q(3));
        fprintf(fid, '%s = %1.16f;\n', '         z', Q(4));
        fprintf(fid, '%s\n', '      }');
        fprintf(fid, '%s\n', '      translation {');
        fprintf(fid, '%s = %1.16f;\n', '         x', x(1));
        fprintf(fid, '%s = %1.16f;\n', '         y', x(2));
        fprintf(fid, '%s = %1.16f;\n', '         z', x(3));
        fprintf(fid, '%s\n', '      }');
    end
    fprintf(fid, '%s\n', '   }');
    fprintf(fid, '%s\n', '}');
    
    fclose(fid);
end
