% ==============================================================
% function ret = findMostFreqDim(list)
% ==============================================================
function ret = findMostFreqDim(list)
    [sortedList Idx]= sort(list);
    countList = zeros(size(sortedList));
    count = 1;
    for i = 1:length(list)
        if i == length(list)
            countList(i) = count;
        elseif sortedList(i) == sortedList(i+1)
            count = count + 1;
        else
            countList(i) = count;
            count = 1;
        end
    end
    [C,id] = max(countList);
    ret = [ list(Idx(id(1))), Idx(id(1))];    
return;
