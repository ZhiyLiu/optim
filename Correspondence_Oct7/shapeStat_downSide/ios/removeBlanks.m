function strVecs = removeBlanks(text)
    [nRow, nCol] = size(text);
    strVecs = [];
    for j = 1:nRow
        line = text(j, :);
        space = isspace(line);
        for k = 1:length(line)
            if space(k) 
                continue;
            else
                break;
            end
        end
        line(1:k-1) = [];
        nLine = zeros(1,nCol);
        for m = 1:length(line)
            nLine(m) = line(m);
        end
        nLine = char(nLine);
        strVecs = [strVecs; nLine];
    end
return;
    
