% Function [M] = mrepLog(Mlist)
%
%   mrepLog(Mlist) computes the logmap of mreps in Mlist.
%
% INPUT : Mlist - cell array of combined(tube, quad) atoms
% OUTPUT : M - cell array of combined(tube, quad) atoms

function M = mrepLog(Mlist)

[row, col] = size(Mlist);

for r = 1:row
    for c = 1:col
    M{r, c} = Log(Mlist{r, c}); 
    end
end

return;